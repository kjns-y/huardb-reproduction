from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import seaborn as sns

from utils import available_gene_symbols, load_config, project_path, setup_logging


def expression_frame(adata: sc.AnnData, genes: list[str]) -> pd.DataFrame:
    selected = available_gene_symbols(adata, genes)
    matrix = adata[:, selected].X
    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    return pd.DataFrame(values, index=adata.obs_names, columns=selected)


def main() -> None:
    parser = argparse.ArgumentParser(description="Explore clone-to-transcriptome, DEG, and V-gene pairing relationships.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    logger = setup_logging("clone_phenotype")
    adata = sc.read_h5ad(project_path("results/objects/cpic_c2_merged.h5ad"))
    clones = pd.read_csv(project_path("results/tables/clonotypes.tsv"), sep="\t")
    hct = pd.read_csv(project_path("results/tables/hct_cells.tsv"), sep="\t")

    min_size = int(config["clonotype"]["transcriptome_correlation_min_size"])
    correlation_ids = clones.loc[clones["clone_size"].ge(min_size), "clonotype_id"].tolist()
    hvg_names = adata.var_names[adata.var["highly_variable"]].tolist()
    clone_means: dict[str, np.ndarray] = {}
    clone_labels = adata.obs["clonotype_id"].astype("string")
    for clone_id in correlation_ids:
        mask = clone_labels.eq(clone_id).fillna(False).to_numpy()
        matrix = adata[mask, hvg_names].X
        clone_means[clone_id] = np.asarray(matrix.mean(axis=0)).ravel()
    if len(clone_means) >= 2:
        mean_frame = pd.DataFrame.from_dict(clone_means, orient="index", columns=hvg_names)
        correlation = mean_frame.T.corr(method="pearson")
    else:
        correlation = pd.DataFrame(index=correlation_ids, columns=correlation_ids, dtype=float)
    correlation.to_csv(project_path("results/tables/clone_transcriptome_correlation.tsv"), sep="\t")
    if not correlation.empty:
        fig, ax = plt.subplots(figsize=(max(6, 0.55 * len(correlation)), max(5, 0.5 * len(correlation))))
        sns.heatmap(correlation.astype(float), vmin=-1, vmax=1, center=0, cmap="vlag", square=True, ax=ax, cbar_kws={"label": "Pearson r"})
        ax.set(title=f"Mean transcriptome correlation among clones with size ≥ {min_size}", xlabel="Clonotype", ylabel="Clonotype")
        fig.tight_layout()
        fig.savefig(project_path("results/figures/clone_transcriptome_correlation.png"), dpi=180, bbox_inches="tight")
        plt.close(fig)

    top_ids = clones.head(10)["clonotype_id"].tolist()
    requested_markers = config["markers"]["clone_expression"]
    marker_values = expression_frame(adata, requested_markers)
    marker_values["clonotype_id"] = clone_labels.to_numpy()
    marker_means = marker_values.loc[marker_values["clonotype_id"].isin(top_ids)].groupby("clonotype_id").mean().reindex(top_ids)
    marker_means.to_csv(project_path("results/tables/top_clone_marker_expression.tsv"), sep="\t")
    fig, ax = plt.subplots(figsize=(max(7, 0.8 * marker_means.shape[1]), 6))
    sns.heatmap(marker_means, cmap="mako", linewidths=0.3, linecolor="white", ax=ax, cbar_kws={"label": "Mean log-normalized expression"})
    ax.set(title="Transcriptional state markers across top clonotypes", xlabel="Gene", ylabel="Clonotype")
    fig.tight_layout()
    fig.savefig(project_path("results/figures/top_clone_marker_expression.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)

    deg_output = project_path("results/tables/clone_deg.tsv")
    if len(clones) >= 2 and clones.iloc[1]["clone_size"] >= 2:
        clone_a, clone_b = clones.iloc[0]["clonotype_id"], clones.iloc[1]["clonotype_id"]
        mask = clone_labels.isin([clone_a, clone_b]).fillna(False).to_numpy()
        subset = adata[mask].copy()
        subset.obs["deg_clone"] = subset.obs["clonotype_id"].astype(str).astype("category")
        sc.tl.rank_genes_groups(subset, groupby="deg_clone", groups=[clone_a], reference=clone_b, method="wilcoxon", pts=True)
        deg = sc.get.rank_genes_groups_df(subset, group=clone_a)
        deg.insert(0, "comparison", f"{clone_a}_vs_{clone_b}")
        deg.to_csv(deg_output, sep="\t", index=False)
        plot = deg.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchanges", "pvals_adj"]).copy()
        plot["minus_log10_padj"] = -np.log10(plot["pvals_adj"].clip(lower=1e-300))
        fig, ax = plt.subplots(figsize=(7, 5))
        significant = plot["pvals_adj"].lt(0.05) & plot["logfoldchanges"].abs().ge(1)
        ax.scatter(plot.loc[~significant, "logfoldchanges"], plot.loc[~significant, "minus_log10_padj"], s=8, c="#B9C1CB", alpha=0.5, linewidths=0)
        ax.scatter(plot.loc[significant, "logfoldchanges"], plot.loc[significant, "minus_log10_padj"], s=10, c="#2A6FBB", alpha=0.7, linewidths=0)
        for _, row in plot.nsmallest(8, "pvals_adj").iterrows():
            ax.text(row["logfoldchanges"], row["minus_log10_padj"], str(row["names"]), fontsize=7)
        ax.axvline(0, color="#333333", linewidth=0.8)
        ax.set(title=f"Exploratory single-cell DEG: {clone_a} vs {clone_b}", xlabel="Log fold change", ylabel="−log10 adjusted p-value")
        fig.tight_layout()
        fig.savefig(project_path("results/figures/clone_deg_volcano.png"), dpi=180, bbox_inches="tight")
        plt.close(fig)
    else:
        pd.DataFrame(columns=["comparison", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"]).to_csv(deg_output, sep="\t", index=False)

    pairing = pd.crosstab(hct["TRA_v"], hct["TRB_v"])
    pairing.to_csv(project_path("results/tables/trav_trbv_pairing.tsv"), sep="\t")
    top_tra = pairing.sum(axis=1).nlargest(20).index
    top_trb = pairing.sum(axis=0).nlargest(20).index
    plot_pairing = pairing.loc[top_tra, top_trb]
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(plot_pairing, cmap="Blues", linewidths=0.2, linecolor="white", ax=ax, cbar_kws={"label": "Strict hcT cells"})
    ax.set(title="Most frequent TRA V × TRB V pairings in strict hcT cells", xlabel="TRB V gene", ylabel="TRA V gene")
    fig.tight_layout()
    fig.savefig(project_path("results/figures/trav_trbv_pairing.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    logger.info("Clone phenotype analysis complete: correlation_clones=%d", len(correlation_ids))


if __name__ == "__main__":
    main()

