from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from utils import load_config, project_path, setup_logging


BLUE = "#2A6FBB"
INK = "#1C2B39"
GREY = "#D8DEE6"
PALETTE = ["#2A6FBB", "#C78B18", "#C66A1B", "#7D8B3A", "#C0527D", "#5C6F82", "#7E57A2", "#3F8C7A", "#A95A4A", "#5865A6"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze clone expansion and map clones onto UMAP.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    load_config(args.config)
    logger = setup_logging("clonotype_analysis")
    adata = sc.read_h5ad(project_path("results/objects/cpic_c2_merged.h5ad"))
    clones = pd.read_csv(project_path("results/tables/clonotypes.tsv"), sep="\t")

    frequencies = clones["clone_size"].to_numpy(dtype=float)
    probabilities = frequencies / frequencies.sum()
    shannon = float(-(probabilities * np.log(probabilities)).sum())
    simpson = float(1 - np.square(probabilities).sum())
    largest = clones.iloc[0]
    tcr_associated = int(adata.obs["has_tcr"].sum())
    strict_hct = int(adata.obs["is_hct"].sum())
    summary = pd.DataFrame([{
        "n_cells": adata.n_obs,
        "tcr_associated_cells": tcr_associated,
        "strict_hct_cells": strict_hct,
        "paired_TRA_TRB_rate": strict_hct / tcr_associated if tcr_associated else np.nan,
        "n_clonotypes": len(clones),
        "n_singletons": int(clones["clone_size"].eq(1).sum()),
        "n_expanded_clones": int(clones["is_expanded"].sum()),
        "largest_clonotype_id": largest["clonotype_id"],
        "largest_clone_size": int(largest["clone_size"]),
        "largest_clone_fraction": float(largest["clone_size"] / strict_hct),
        "shannon_diversity": shannon,
        "simpson_diversity": simpson,
        "hill_q0": len(clones),
        "hill_q1": float(np.exp(shannon)),
        "hill_q2": float(1 / np.square(probabilities).sum()),
    }])
    summary.to_csv(project_path("results/tables/analysis_summary.tsv"), sep="\t", index=False)
    clones.groupby("clone_size").size().rename("n_clonotypes").reset_index().to_csv(
        project_path("results/tables/clone_size_distribution.tsv"), sep="\t", index=False
    )

    top = clones.head(10).sort_values("clone_size", ascending=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(top["clonotype_id"], top["clone_size"], color=BLUE, edgecolor=INK)
    ax.set(xlabel="Cells in clonotype", ylabel="Clonotype", title="Top 10 CPIc_C2 clonotypes by exact paired CDR3 nucleotide sequence")
    for index, value in enumerate(top["clone_size"]):
        ax.text(value, index, f" {value}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(project_path("results/figures/top10_clonotypes.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)

    coords = adata.obsm["X_umap"]
    clone_ids = adata.obs["clonotype_id"].astype("string")
    top_ids = clones.head(10)["clonotype_id"].tolist()
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(coords[:, 0], coords[:, 1], s=7, c=GREY, alpha=0.45, linewidths=0)
    for color, clone_id in zip(PALETTE, top_ids):
        mask = clone_ids.eq(clone_id).fillna(False).to_numpy()
        ax.scatter(coords[mask, 0], coords[mask, 1], s=20, c=color, label=f"{clone_id} (n={mask.sum()})", edgecolors=INK, linewidths=0.2)
    ax.set(title="Top 10 clonotypes projected onto CPIc_C2 transcriptome", xlabel="UMAP1", ylabel="UMAP2")
    ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(project_path("results/figures/umap_top_clones.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)

    top_one = top_ids[0]
    mask = clone_ids.eq(top_one).fillna(False).to_numpy()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(coords[:, 0], coords[:, 1], s=7, c=GREY, alpha=0.45, linewidths=0)
    ax.scatter(coords[mask, 0], coords[mask, 1], s=24, c=BLUE, edgecolors=INK, linewidths=0.25)
    ax.set(title=f"{top_one} projected onto CPIc_C2 transcriptome (n={mask.sum()})", xlabel="UMAP1", ylabel="UMAP2")
    fig.tight_layout()
    fig.savefig(project_path("results/figures/umap_clone_01.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)

    hct_obs = adata.obs.loc[adata.obs["is_hct"]].copy()
    contingency = pd.crosstab(hct_obs["clonotype_id"], hct_obs["leiden"])
    contingency.to_csv(project_path("results/tables/clone_by_cluster.tsv"), sep="\t")
    heat = contingency.reindex(top_ids).fillna(0)
    row_sum = heat.sum(axis=1).replace(0, np.nan)
    heat_prop = heat.div(row_sum, axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(max(7, 0.55 * heat.shape[1]), 6))
    sns.heatmap(heat_prop, cmap="Blues", linewidths=0.3, linecolor="white", ax=ax, cbar_kws={"label": "Within-clone fraction"})
    ax.set(title="Top clonotype distribution across Leiden clusters", xlabel="Leiden cluster", ylabel="Clonotype")
    fig.tight_layout()
    fig.savefig(project_path("results/figures/clone_by_cluster.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    logger.info("Clone analysis: clones=%d expanded=%d largest=%d", len(clones), int(clones["is_expanded"].sum()), int(largest["clone_size"]))


if __name__ == "__main__":
    main()

