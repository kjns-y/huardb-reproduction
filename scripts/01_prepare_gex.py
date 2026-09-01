from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy import sparse

from utils import load_config, project_path, setup_logging


def main() -> None:
    parser = argparse.ArgumentParser(description="Read the CPIc_C2 processed 10x GEX matrix.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    logger = setup_logging("prepare_gex")

    source = project_path(config["data"]["gex_source_dir"])
    output = project_path("results/objects/cpic_c2_raw.h5ad")
    output.parent.mkdir(parents=True, exist_ok=True)
    if not source.exists():
        raise FileNotFoundError(f"GEX source directory does not exist: {source}")

    adata = sc.read_10x_mtx(source, var_names="gene_symbols", make_unique=False, cache=False)
    raw_symbols = pd.Index(adata.var_names.astype(str))
    duplicate_genes = int(raw_symbols.duplicated().sum())
    duplicate_barcodes = int(adata.obs_names.duplicated().sum())
    adata.var["gene_symbol"] = raw_symbols.to_numpy()
    adata.var_names_make_unique(join="-dup")
    adata.obs_names = adata.obs_names.astype(str)
    adata.uns["source_audit"] = {
        "sample_id": config["project"]["sample_id"],
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "duplicate_barcodes": duplicate_barcodes,
        "duplicate_gene_symbols": duplicate_genes,
        "is_sparse": bool(sparse.issparse(adata.X)),
        "matrix_orientation": "cells_x_genes",
    }
    if duplicate_barcodes:
        raise ValueError(f"GEX contains {duplicate_barcodes} duplicate barcodes")
    if not sparse.issparse(adata.X):
        raise ValueError("GEX matrix unexpectedly became dense")
    adata.write_h5ad(output, compression="gzip")
    logger.info("GEX loaded: cells=%d genes=%d sparse=%s", adata.n_obs, adata.n_vars, True)


if __name__ == "__main__":
    main()

