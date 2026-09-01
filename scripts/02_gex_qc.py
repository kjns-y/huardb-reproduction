from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from utils import load_config, project_path, setup_logging


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate and apply GEX QC for CPIc_C2.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    logger = setup_logging("gex_qc")

    adata = sc.read_h5ad(project_path("results/objects/cpic_c2_raw.h5ad"))
    gene_symbols = adata.var["gene_symbol"].astype(str)
    adata.var["mt"] = gene_symbols.str.upper().str.startswith("MT-").to_numpy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    initial = pd.DataFrame([{
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "median_counts": float(np.median(adata.obs["total_counts"])),
        "median_genes": float(np.median(adata.obs["n_genes_by_counts"])),
        "duplicate_barcodes": int(adata.uns["source_audit"]["duplicate_barcodes"]),
        "duplicate_genes": int(adata.uns["source_audit"]["duplicate_gene_symbols"]),
        "sparse_matrix": bool(sparse.issparse(adata.X)),
    }])
    initial.to_csv(project_path("results/qc/gex_initial_qc.tsv"), sep="\t", index=False)

    figures = project_path("results/figures")
    figures.mkdir(parents=True, exist_ok=True)
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.25,
        multi_panel=True,
        show=False,
    )
    plt.gcf().suptitle("CPIc_C2 GEX QC before filtering", y=1.02)
    plt.savefig(figures / "gex_qc_violin.png", dpi=180, bbox_inches="tight")
    plt.close("all")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"], s=5, alpha=0.35, color="#2A6FBB")
    axes[0].set(xlabel="Total counts", ylabel="Genes detected", title="Library size vs genes")
    axes[1].scatter(adata.obs["total_counts"], adata.obs["pct_counts_mt"], s=5, alpha=0.35, color="#C66A1B")
    axes[1].axhline(config["gex_qc"]["max_pct_mt"], color="#333333", linestyle="--", linewidth=1)
    axes[1].set(xlabel="Total counts", ylabel="Mitochondrial counts (%)", title="Library size vs mitochondrial fraction")
    fig.suptitle("CPIc_C2 GEX QC before filtering")
    fig.tight_layout()
    fig.savefig(figures / "gex_qc_scatter.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    steps: list[dict[str, int | str]] = []
    before = adata.n_obs
    keep_genes = adata.obs["n_genes_by_counts"].ge(config["gex_qc"]["min_genes"])
    after = int(keep_genes.sum())
    steps.append({"step": "n_genes_by_counts>=200", "cells_before": before, "cells_after": after, "cells_removed": before - after})
    adata = adata[keep_genes].copy()
    before = adata.n_obs
    keep_mt = adata.obs["pct_counts_mt"].le(config["gex_qc"]["max_pct_mt"])
    after = int(keep_mt.sum())
    steps.append({"step": "pct_counts_mt<=20", "cells_before": before, "cells_after": after, "cells_removed": before - after})
    adata = adata[keep_mt].copy()
    steps.append({"step": "doublet_removal_not_run", "cells_before": adata.n_obs, "cells_after": adata.n_obs, "cells_removed": 0})
    pd.DataFrame(steps).to_csv(project_path("results/qc/gex_filtering_audit.tsv"), sep="\t", index=False)
    adata.uns["gex_qc"] = {
        "min_genes": int(config["gex_qc"]["min_genes"]),
        "max_pct_mt": float(config["gex_qc"]["max_pct_mt"]),
        "doublet_method": config["gex_qc"]["doublet_method"],
    }
    adata.write_h5ad(project_path("results/objects/cpic_c2_qc.h5ad"), compression="gzip")
    logger.info("GEX QC: before=%d after=%d", int(initial.loc[0, "n_cells"]), adata.n_obs)


if __name__ == "__main__":
    main()

