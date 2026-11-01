#!/usr/bin/env python3
"""Stream a Seurat v5 .rds into an .h5ad without R, and without ever holding
the whole object in memory.

Why this exists: sc_UMI_count.rds expands to ~5 GB, the Cowork VM has 4 GB of
RAM and no R, so readRDS is not an option at all. This walks R's serialization
format directly and pipes the two large sparse-matrix arrays straight into
HDF5, so peak memory stays at one chunk.

    python3 rds_to_h5ad.py scan    in.rds  manifest.json
    python3 rds_to_h5ad.py extract in.rds  manifest.json  out.h5ad

Two passes on purpose. The scan records where the big vectors are and what
tag each sits under; the extract re-walks the same file deterministically and
streams only the targeted ones. That avoids threading context through the
parser and keeps the scan cheap enough to re-run.

A dgCMatrix is genes x cells in CSC. CSC(genes x cells) and CSR(cells x genes)
are the same three arrays, so the transpose is free: i/p/x are written
unchanged as a CSR cells x genes matrix, which is the orientation AnnData
wants.
"""
from __future__ import annotations

import gzip
import json
import struct
import sys

import numpy as np

BIG = 200_000          # elements; above this a vector is streamed, not kept
CHUNK = 1 << 22        # 4 MiB

# R serialization type codes (src/main/serialize.c)
NILSXP, SYMSXP, LISTSXP, CLOSXP, ENVSXP, PROMSXP, LANGSXP = 0, 1, 2, 3, 4, 5, 6
CHARSXP, LGLSXP, INTSXP, REALSXP, CPLXSXP, STRSXP = 9, 10, 13, 14, 15, 16
DOTSXP, VECSXP, EXPRSXP, BCODESXP, EXTPTRSXP = 17, 19, 20, 21, 22
WEAKREFSXP, RAWSXP, S4SXP = 23, 24, 25
ALTREP_SXP, ATTRLISTSXP, ATTRLANGSXP = 238, 239, 240
BASEENV_SXP, EMPTYENV_SXP, MISSINGARG_SXP = 241, 242, 251
GENERICREFSXP, CLASSREFSXP, PERSISTSXP, PACKAGESXP, NAMESPACESXP = 245, 246, 247, 248, 249
BASENAMESPACE_SXP, UNBOUNDVALUE_SXP, GLOBALENV_SXP, NILVALUE_SXP, REFSXP = 250, 252, 253, 254, 255
BCREPREF, BCREPDEF = 243, 244

PAIRLIKE = {LISTSXP, LANGSXP, CLOSXP, PROMSXP, DOTSXP}


class Sym(str):
    """An R symbol. Subclasses str so tag names compare naturally."""


class Env:
    __slots__ = ("frame", "hashtab", "enclos", "attr")


class Big:
    """Placeholder for a vector too large to keep. `idx` is its ordinal in the
    walk, which is what the extract pass matches on."""
    __slots__ = ("idx", "rtype", "n", "pos", "path")

    def __init__(self, idx, rtype, n, pos, path):
        self.idx, self.rtype, self.n, self.pos, self.path = idx, rtype, n, pos, path

    def __repr__(self):
        return f"<Big #{self.idx} {self.rtype} n={self.n}>"


class Reader:
    def __init__(self, path):
        self.f = gzip.open(path, "rb")
        self.pos = 0
        self.refs = []
        self.bigs = []
        self.targets = {}      # idx -> callable(reader, dtype, itemsize, n)
        self.n_big = 0

    # --- primitives ---
    def read(self, n):
        out = self.f.read(n)
        if len(out) != n:
            raise EOFError(f"wanted {n}, got {len(out)} at {self.pos}")
        self.pos += n
        return out

    def skip(self, n):
        left = n
        while left > 0:
            self.read(min(left, CHUNK))
            left -= min(left, CHUNK)

    def i32(self):
        return struct.unpack(">i", self.read(4))[0]

    def vlen(self):
        n = self.i32()
        if n == -1:                       # long vector
            hi, lo = self.i32(), self.i32()
            return (hi << 32) | (lo & 0xFFFFFFFF)
        return n


def read_header(r):
    magic = r.read(2)
    if magic != b"X\n":
        raise ValueError(f"not an XDR-format RDS (magic {magic!r}); "
                         "a version-2 or ascii .rds is not handled here")
    ver, writer, minreader = r.i32(), r.i32(), r.i32()
    enc = None
    if ver >= 3:
        enc = r.read(r.i32()).decode("ascii")
    return dict(version=ver, writer=writer, min_reader=minreader, encoding=enc)


def atomic(r, rtype, dtype, itemsize, n, path):
    """Either keep a small vector, stream a targeted big one, or skip it."""
    if n <= BIG:
        return np.frombuffer(r.read(n * itemsize), dtype=dtype).astype(
            np.int32 if itemsize == 4 else np.float64)
    idx = r.n_big
    r.n_big += 1
    b = Big(idx, rtype, n, r.pos, path)
    r.bigs.append(b)
    sink = r.targets.get(idx)
    if sink is None:
        r.skip(n * itemsize)
    else:
        sink(r, dtype, itemsize, n)
    return b


def read_bc(r, path):
    """Seurat stores SeuratCommand objects whose params include functions, and
    a byte-compiled function serializes as BCODESXP. We only need to advance
    the stream past it correctly, so the result is discarded -- but the
    traversal has to match R's ReadBC exactly or every later offset is wrong."""
    nreps = r.i32()
    return read_bc1(r, [None] * (nreps + 1), path)


def read_bc1(r, reps, path):
    read_item(r, path + ("<bc-code>",))          # the code, an INTSXP
    read_bc_consts(r, reps, path)
    return {"__bytecode__": True}


def read_bc_consts(r, reps, path):
    n = r.i32()
    for _ in range(n):
        t = r.i32()
        if t == BCODESXP:
            read_bc1(r, reps, path)
        elif t in (LANGSXP, LISTSXP, BCREPDEF, BCREPREF, ATTRLANGSXP, ATTRLISTSXP):
            read_bc_lang(r, t, reps, path)
        else:
            # the type int was a pad; the item reads its own flags
            read_item(r, path + ("<bc-const>",))


def read_bc_lang(r, t, reps, path):
    if t == BCREPREF:
        return reps[r.i32()]
    if t in (BCREPDEF, LANGSXP, LISTSXP, ATTRLANGSXP, ATTRLISTSXP):
        pos = -1
        if t == BCREPDEF:
            pos = r.i32()
            t = r.i32()
        has_attr = False
        if t == ATTRLANGSXP:
            t, has_attr = LANGSXP, True
        elif t == ATTRLISTSXP:
            t, has_attr = LISTSXP, True
        node = {"__lang__": True}
        if pos >= 0:
            reps[pos] = node
        if has_attr:
            read_item(r, path + ("<bc-attr>",))
        read_item(r, path + ("<bc-tag>",))
        read_bc_lang(r, r.i32(), reps, path)
        read_bc_lang(r, r.i32(), reps, path)
        return node
    return read_item(r, path + ("<bc-lang>",))


def read_item(r, path=()):
    flags = r.i32()
    t = flags & 0xFF
    has_attr = bool((flags >> 9) & 1)
    has_tag = bool((flags >> 10) & 1)

    if t == ATTRLISTSXP:
        t, has_attr = LISTSXP, True
    elif t == ATTRLANGSXP:
        t, has_attr = LANGSXP, True

    # --- singletons and references ---
    if t == NILVALUE_SXP or t == NILSXP:
        return None
    if t == REFSXP:
        i = flags >> 8
        if i == 0:
            i = r.i32()
        return r.refs[i - 1]
    if t in (GLOBALENV_SXP, EMPTYENV_SXP, BASEENV_SXP, BASENAMESPACE_SXP,
             MISSINGARG_SXP, UNBOUNDVALUE_SXP):
        return Sym(f".{t}")

    if t == SYMSXP:
        name = read_item(r, path)
        s = Sym(name if isinstance(name, str) else "")
        r.refs.append(s)
        return s

    if t in (PACKAGESXP, NAMESPACESXP):
        val = read_item(r, path)
        r.refs.append(val)
        return val

    if t == ENVSXP:
        e = Env()
        r.refs.append(e)              # must register before reading contents
        r.i32()                       # locked
        e.enclos = read_item(r, path + ("<enclos>",))
        e.frame = read_item(r, path + ("<frame>",))
        e.hashtab = read_item(r, path + ("<hash>",))
        e.attr = read_item(r, path + ("<attr>",))
        return e

    if t in PAIRLIKE:
        attr = read_item(r, path + ("<attr>",)) if has_attr else None
        tag = read_item(r, path) if has_tag else None
        name = str(tag) if tag is not None else "?"
        car = read_item(r, path + (name,))
        cdr = read_item(r, path)
        node = {"__tag__": name, "__car__": car, "__cdr__": cdr, "__attr__": attr}
        return node

    if t == CHARSXP:
        n = r.i32()
        if n == -1:
            return None
        return r.read(n).decode("utf-8", "replace")

    # --- vectors ---
    if t in (LGLSXP, INTSXP):
        val = atomic(r, "int", ">i4", 4, r.vlen(), path)
    elif t == REALSXP:
        val = atomic(r, "real", ">f8", 8, r.vlen(), path)
    elif t == CPLXSXP:
        n = r.vlen(); r.skip(n * 16); val = f"<complex n={n}>"
    elif t == RAWSXP:
        n = r.vlen(); r.skip(n); val = f"<raw n={n}>"
    elif t == STRSXP:
        n = r.vlen()
        if n > BIG:
            idx = r.n_big; r.n_big += 1
            b = Big(idx, "str", n, r.pos, path); r.bigs.append(b)
            for _ in range(n):
                read_item(r, path)
            val = b
        else:
            val = [read_item(r, path) for _ in range(n)]
    elif t in (VECSXP, EXPRSXP):
        n = r.vlen()
        val = [read_item(r, path + (f"[{i}]",)) for i in range(n)]
    elif t == S4SXP:
        val = {"__s4__": True}
    elif t == ALTREP_SXP:
        info = read_item(r, path + ("<altrep-info>",))
        state = read_item(r, path + ("<altrep-state>",))
        read_item(r, path + ("<altrep-attr>",))
        val = {"__altrep__": True, "info": info, "state": state}
    elif t == BCODESXP:
        val = read_bc(r, path)
    elif t == EXTPTRSXP:
        read_item(r, path + ("<extptr-prot>",))
        read_item(r, path + ("<extptr-tag>",))
        val = "<externalptr>"
    elif t == WEAKREFSXP:
        val = "<weakref>"
    elif t in (CLASSREFSXP, GENERICREFSXP, PERSISTSXP):
        val = read_item(r, path + ("<ref>",))
    else:
        raise NotImplementedError(f"unhandled SEXP type {t} at {'/'.join(path)}")

    if has_attr:
        a = read_item(r, path + ("<attr>",))
        val = {"__val__": val, "__attr__": a}
    return val


# --- helpers for poking at the parsed tree -------------------------------
def plist(node):
    """Flatten an R pairlist into {tag: value}."""
    out = {}
    while isinstance(node, dict) and "__tag__" in node:
        out[node["__tag__"]] = node["__car__"]
        node = node["__cdr__"]
    return out


def attrs_of(x):
    if isinstance(x, dict) and "__attr__" in x and "__val__" in x:
        return plist(x["__attr__"])
    if isinstance(x, dict) and "__attr__" in x:
        return plist(x["__attr__"])
    return {}


def value_of(x):
    return x["__val__"] if isinstance(x, dict) and "__val__" in x else x


def class_of(x):
    c = attrs_of(x).get("class")
    c = value_of(c)
    if isinstance(c, list) and c:
        return str(c[0])
    return str(c) if isinstance(c, str) else None


def walk(x, depth=0, path="", out=None, maxdepth=7):
    """Yield (path, class, kind) for every node, for eyeballing the layout."""
    if out is None:
        out = []
    if depth > maxdepth:
        return out
    cls = class_of(x)
    v = value_of(x)
    kind = type(v).__name__
    if isinstance(v, Big):
        kind = f"Big(n={v.n},{v.rtype})"
    out.append((path, cls, kind))
    a = attrs_of(x)
    for k, sub in a.items():
        if k == "class":
            continue
        walk(sub, depth + 1, f"{path}@{k}", out, maxdepth)
    if isinstance(v, list):
        names = value_of(a.get("names")) or []
        for i, sub in enumerate(v):
            nm = names[i] if i < len(names) and isinstance(names[i], str) else f"[{i}]"
            walk(sub, depth + 1, f"{path}${nm}", out, maxdepth)
    elif isinstance(v, dict) and "__tag__" in v:
        for k, sub in plist(v).items():
            walk(sub, depth + 1, f"{path}${k}", out, maxdepth)
    return out


def find_dgc(root):
    """Locate every dgCMatrix in the tree, with its slots."""
    found = []

    def rec(x, path, depth=0):
        if depth > 9:
            return
        a = attrs_of(x)
        if class_of(x) in ("dgCMatrix", "dgTMatrix", "dgRMatrix"):
            found.append((path, class_of(x), a))
        v = value_of(x)
        for k, sub in a.items():
            if k != "class":
                rec(sub, f"{path}@{k}", depth + 1)
        if isinstance(v, list):
            names = value_of(a.get("names")) or []
            for i, sub in enumerate(v):
                nm = names[i] if i < len(names) and isinstance(names[i], str) else f"[{i}]"
                rec(sub, f"{path}${nm}", depth + 1)
        elif isinstance(v, dict) and "__tag__" in v:
            for k, sub in plist(v).items():
                rec(sub, f"{path}${k}", depth + 1)

    rec(root, "")
    return found


# ------------------------------------------------------------------ scan
def cmd_scan(rds, manifest_path):
    r = Reader(rds)
    hdr = read_header(r)
    print("header:", hdr)
    root = read_item(r)
    print(f"walked {r.pos / 1e9:.2f} GB; {len(r.bigs)} big vectors\n")

    print("=== top-level layout ===")
    for p, cls, kind in walk(root, maxdepth=4):
        if p:
            print(f"  {p:<58} {str(cls):<14} {kind}")

    mats = find_dgc(root)
    print(f"\n=== sparse matrices found: {len(mats)} ===")
    out = {"header": hdr, "matrices": []}
    for path, cls, a in mats:
        dim = value_of(a.get("Dim"))
        dim = [int(x) for x in np.asarray(dim).ravel()] if dim is not None else None
        rec = {"path": path, "class": cls, "dim": dim, "slots": {}}
        for slot in ("i", "p", "x"):
            s = value_of(a.get(slot))
            if isinstance(s, Big):
                rec["slots"][slot] = {"big_idx": s.idx, "n": s.n, "rtype": s.rtype}
            elif s is not None:
                arr = np.asarray(s).ravel()
                rec["slots"][slot] = {"inline": arr.tolist(), "n": int(arr.size)}
        dn = value_of(a.get("Dimnames"))
        if isinstance(dn, list):
            rec["dimnames"] = []
            for d in dn:
                dv = value_of(d)
                if isinstance(dv, Big):
                    rec["dimnames"].append({"big_idx": dv.idx, "n": dv.n})
                elif isinstance(dv, list):
                    rec["dimnames"].append({"inline": [str(z) for z in dv], "n": len(dv)})
                else:
                    rec["dimnames"].append(None)
        out["matrices"].append(rec)
        print(f"  {path}  {cls}  dim={dim}")
        for k, v in rec["slots"].items():
            print(f"      {k}: {v if 'inline' not in v else {'n': v['n']}}")
        for j, d in enumerate(rec.get("dimnames") or []):
            print(f"      Dimnames[{j}]: {d if d is None or 'inline' not in d else {'n': d['n']}}")

    out["bigs"] = [{"idx": b.idx, "rtype": b.rtype, "n": b.n,
                    "path": "/".join(b.path)[-120:]} for b in r.bigs]
    print(f"\n=== the 8 largest streamed vectors ===")
    for b in sorted(r.bigs, key=lambda z: -z.n)[:8]:
        print(f"  #{b.idx:<3} {b.rtype:<5} n={b.n:>12,}  {'/'.join(b.path)[-90:]}")

    with open(manifest_path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {manifest_path}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "scan"
    if cmd == "scan":
        cmd_scan(sys.argv[2], sys.argv[3])
    else:
        from rds_extract import cmd_extract      # kept separate, see scripts/
        cmd_extract(*sys.argv[2:])
