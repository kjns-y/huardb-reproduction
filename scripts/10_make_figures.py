from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from utils import available_gene_symbols, load_config, project_path, setup_logging


def main() -> None:
    parser = argparse.ArgumentParser(description="Create marker UMAP panels and verify core figure inventory.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    logger = setup_logging("make_figures")
    adata = sc.read_h5ad(project_path("results/objects/cpic_c2_merged.h5ad"))
    markers = available_gene_symbols(adata, config["markers"]["feature_plot"])
    missing_markers = sorted(set(config["markers"]["feature_plot"]) - set(markers))
    if not markers:
        raise ValueError("No requested marker genes found in GEX matrix")
    sc.pl.umap(
        adata,
        color=markers,
        ncols=3,
        cmap="viridis",
        frameon=True,
        wspace=0.35,
        show=False,
    )
    plt.gcf().suptitle("CPIc_C2 canonical T-cell marker expression", y=1.01)
    plt.savefig(project_path("results/figures/umap_marker_genes.png"), dpi=180, bbox_inches="tight")
    plt.close("all")

    expected = [
        "gex_qc_violin.png",
        "gex_qc_scatter.png",
        "tcr_chain_configuration.png",
        "umap_leiden.png",
        "umap_marker_genes.png",
        "top10_clonotypes.png",
        "umap_top_clones.png",
        "umap_clone_01.png",
        "clone_by_cluster.png",
        "top_clone_marker_expression.png",
        "trav_trbv_pairing.png",
    ]
    rows = []
    for name in expected:
        path = project_path("results/figures") / name
        rows.append({"figure": name, "exists": path.exists(), "size_bytes": path.stat().st_size if path.exists() else 0})
    pd.DataFrame(rows).to_csv(project_path("results/qc/figure_inventory.tsv"), sep="\t", index=False)
    pd.DataFrame({"missing_requested_marker": missing_markers}).to_csv(
        project_path("results/qc/missing_marker_genes.tsv"), sep="\t", index=False
    )
    if not all(row["exists"] and row["size_bytes"] > 0 for row in rows):
        raise RuntimeError("One or more required figures are missing or empty")
    logger.info("Core figures verified: %d; missing requested markers=%d", len(rows), len(missing_markers))


if __name__ == "__main__":
    main()

