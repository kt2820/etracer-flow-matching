#!/usr/bin/env python3
"""Second pass: stream the counts dgCMatrix out of the Seurat .rds into an
.h5ad, using the manifest written by `rds_to_h5ad.py scan`.

    python3 rds_extract.py in.rds manifest.json out.h5ad

Peak memory is one chunk (~32 MB), so this runs in a 4 GB VM against a 5 GB
object. Nothing is held whole.

A dgCMatrix is CSC over (genes x cells): `i` holds gene row indices and `p`
holds per-cell offsets. Those are byte-for-byte the CSR arrays of the
transpose, (cells x genes), which is the orientation AnnData wants -- so the
transpose costs nothing and no reindexing happens here.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import h5py

_spec = importlib.util.spec_from_file_location(
    "rdsmod", str(Path(__file__).resolve().parent / "rds_to_h5ad.py"))
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)

STEP = 1 << 21          # elements per chunk (2 M -> 16 MB of float64 in flight)


def find_string_vec(root, want_path, want_n):
    """Pull a names vector out of the parsed tree, by path if possible and by
    length as a fallback."""
    best = None
    def rec(x, path, depth=0):
        nonlocal best
        if depth > 10 or best is not None:
            return
        a, v = R.attrs_of(x), R.value_of(x)
        if (isinstance(v, list) and len(v) == want_n and v
                and isinstance(v[0], str) and path == want_path):
            best = [str(z) for z in v]
            return
        for k, sub in a.items():
            if k != "class":
                rec(sub, f"{path}@{k}", depth + 1)
        if isinstance(v, list):
            names = R.value_of(a.get("names")) or []
            for i, sub in enumerate(v):
                nm = (names[i] if i < len(names) and isinstance(names[i], str)
                      else f"[{i}]")
                rec(sub, f"{path}${nm}", depth + 1)
        elif isinstance(v, dict) and "__tag__" in v:
            for k, sub in R.plist(v).items():
                rec(sub, f"{path}${k}", depth + 1)
    rec(root, "")
    if best is None:                       # fall back: first vector of that length
        hits = []
        def rec2(x, path, depth=0):
            a, v = R.attrs_of(x), R.value_of(x)
            if depth > 10 or hits:
                return
            if isinstance(v, list) and len(v) == want_n and v and isinstance(v[0], str):
                hits.append([str(z) for z in v]); return
            for k, sub in a.items():
                if k != "class":
                    rec2(sub, f"{path}@{k}", depth + 1)
            if isinstance(v, list):
                for i, sub in enumerate(v):
                    rec2(sub, f"{path}[{i}]", depth + 1)
            elif isinstance(v, dict) and "__tag__" in v:
                for k, sub in R.plist(v).items():
                    rec2(sub, f"{path}${k}", depth + 1)
        rec2(root, "")
        best = hits[0] if hits else None
    return best


def make_sink(dset, out_dtype, label, check=None):
    t0 = time.time()
    state = {"max": -np.inf, "nonint": 0}

    def sink(r, dtype, itemsize, n):
        done = 0
        while done < n:
            k = min(STEP, n - done)
            raw = np.frombuffer(r.read(k * itemsize), dtype=dtype)
            arr = raw.astype(out_dtype)
            if check == "int":
                if arr.size:
                    state["max"] = max(state["max"], float(raw.max()))
                state["nonint"] += int(np.count_nonzero(raw != np.floor(raw)))
            elif check == "idx" and arr.size:
                state["max"] = max(state["max"], float(raw.max()))
            dset[done:done + k] = arr
            done += k
            if done % (50 * STEP) < STEP or done == n:
                el = time.time() - t0
                print(f"    {label}: {done:>12,}/{n:,} "
                      f"({done/n:5.1%})  {done*itemsize/1e6/max(el,.01):6.1f} MB/s in",
                      flush=True)
        state["n"] = n
    return sink, state


def cmd_extract(rds, manifest_path, out_h5ad, which=0):
    man = json.load(open(manifest_path))
    mat = man["matrices"][int(which)]
    n_genes, n_cells = mat["dim"]
    nnz = mat["slots"]["i"]["n"]
    p = np.asarray(mat["slots"]["p"]["inline"], dtype=np.int64)
    i_idx = mat["slots"]["i"]["big_idx"]
    x_idx = mat["slots"]["x"]["big_idx"]

    print(f"{mat['path']}  {mat['class']}")
    print(f"  {n_genes:,} genes x {n_cells:,} cells, nnz = {nnz:,} "
          f"({nnz/(n_genes*n_cells):.1%} dense)")
    print(f"  CSC(genes x cells) reused as CSR(cells x genes): no transpose\n")
    assert p[0] == 0 and p[-1] == nnz and len(p) == n_cells + 1
    assert nnz < 2**31 - 1, "indptr would overflow int32"

    out = Path(out_h5ad)
    f = h5py.File(out, "w")
    f.attrs["encoding-type"] = "anndata"
    f.attrs["encoding-version"] = "0.1.0"

    X = f.create_group("X")
    X.attrs["encoding-type"] = "csr_matrix"
    X.attrs["encoding-version"] = "0.1.0"
    X.attrs["shape"] = np.array([n_cells, n_genes], dtype="int64")
    X.create_dataset("indptr", data=p.astype(np.int32), compression="gzip",
                     compression_opts=1)
    d_data = X.create_dataset("data", shape=(nnz,), dtype="float32",
                              chunks=(min(1 << 20, nnz),),
                              compression="gzip", compression_opts=1)
    d_ind = X.create_dataset("indices", shape=(nnz,), dtype="int32",
                             chunks=(min(1 << 20, nnz),),
                             compression="gzip", compression_opts=1)

    sink_i, st_i = make_sink(d_ind, np.int32, "indices (gene ids)", check="idx")
    sink_x, st_x = make_sink(d_data, np.float32, "data (UMI counts)", check="int")

    r = R.Reader(rds)
    r.targets = {i_idx: sink_i, x_idx: sink_x}
    hdr = R.read_header(r)
    t0 = time.time()
    root = R.read_item(r)
    print(f"\n  walked {r.pos/1e9:.2f} GB in {time.time()-t0:.0f} s")

    print(f"  max gene index {st_i['max']:.0f} (must be < {n_genes})")
    print(f"  max UMI value  {st_x['max']:.0f};  non-integer entries: {st_x['nonint']}")
    assert st_i["max"] < n_genes, "gene index out of range"
    assert st_x["nonint"] == 0, "counts are not integers -- not a raw count layer"

    genes = find_string_vec(root, "@assays$RNA@features@dimnames$[0]", n_genes)
    cells = find_string_vec(root, "@assays$RNA@cells@dimnames$[0]", n_cells)
    assert genes and len(genes) == n_genes, "gene names not recovered"
    assert cells and len(cells) == n_cells, "cell names not recovered"
    print(f"  names: {len(cells):,} cells (e.g. {cells[0]}), "
          f"{len(genes):,} genes (e.g. {genes[0]})")

    sdt = h5py.string_dtype(encoding="utf-8")
    for name, vals in (("obs", cells), ("var", genes)):
        g = f.create_group(name)
        g.attrs["encoding-type"] = "dataframe"
        g.attrs["encoding-version"] = "0.2.0"
        g.attrs["_index"] = "_index"
        g.attrs["column-order"] = np.array([], dtype=sdt)
        ds = g.create_dataset("_index", data=np.array(vals, dtype=object),
                              dtype=sdt, compression="gzip", compression_opts=1)
        ds.attrs["encoding-type"] = "string-array"
        ds.attrs["encoding-version"] = "0.2.0"

    for name in ("layers", "obsm", "varm", "obsp", "varp", "uns"):
        g = f.create_group(name)
        g.attrs["encoding-type"] = "dict"
        g.attrs["encoding-version"] = "0.1.0"

    f.close()
    mb = out.stat().st_size / 1e6
    print(f"\nwrote {out}  ({mb:,.0f} MB on disk, "
          f"{nnz*8/1e6:,.0f} MB of raw float32+int32 before compression)")


if __name__ == "__main__":
    cmd_extract(*sys.argv[1:])
