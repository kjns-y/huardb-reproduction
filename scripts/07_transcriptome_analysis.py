from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from utils import available_gene_symbols, load_config, project_path, setup_logging


def main() -> None:
    parser = argparse.ArgumentParser(description="Run deterministic Scanpy preprocessing and exploratory T-cell state annotation.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    seed = int(config["project"]["random_seed"])
    logger = setup_logging("transcriptome")
    np.random.seed(seed)

    adata = sc.read_h5ad(project_path("results/objects/cpic_c2_merged_preprocessed.h5ad"))
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=float(config["transcriptome"]["target_sum"]))
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(int(config["transcriptome"]["n_top_genes"]), adata.n_vars),
        flavor="seurat",
        inplace=True,
    )
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    n_comps = min(int(config["transcriptome"]["n_pcs"]), hvg.n_obs - 1, hvg.n_vars - 1)
    sc.tl.pca(hvg, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"].copy()
    adata.uns["pca_variance"] = hvg.uns["pca"]["variance"]
    adata.uns["pca_variance_ratio"] = hvg.uns["pca"]["variance_ratio"]
    sc.pp.neighbors(
        adata,
        n_neighbors=int(config["transcriptome"]["n_neighbors"]),
        n_pcs=n_comps,
        random_state=seed,
    )
    sc.tl.leiden(
        adata,
        resolution=float(config["transcriptome"]["leiden_resolution"]),
        random_state=seed,
        key_added="leiden",
    )
    sc.tl.umap(adata, random_state=seed)

    score_columns: list[str] = []
    for state, requested in config["markers"]["state_signatures"].items():
        genes = available_gene_symbols(adata, requested)
        if not genes:
            logger.warning("No genes available for state signature %s", state)
            continue
        column = f"state_score_{state}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=column, random_state=seed, use_raw=False)
        score_columns.append(column)
    if score_columns:
        scores = adata.obs[score_columns]
        adata.obs["marker_state"] = scores.idxmax(axis=1).str.replace("state_score_", "", regex=False).astype("category")
    else:
        adata.obs["marker_state"] = pd.Categorical(["unassigned"] * adata.n_obs)

    adata.uns["transcriptome_method"] = {
        "normalization": "normalize_total + log1p",
        "hvg": int(adata.var["highly_variable"].sum()),
        "scaled_for_pca_only": True,
        "n_pcs": int(n_comps),
        "n_neighbors": int(config["transcriptome"]["n_neighbors"]),
        "leiden_resolution": float(config["transcriptome"]["leiden_resolution"]),
        "random_seed": seed,
        "phenotype_annotation": "exploratory multi-gene marker scores",
    }
    output = project_path("results/objects/cpic_c2_merged.h5ad")
    adata.write_h5ad(output, compression="gzip")

    sc.pl.umap(adata, color="leiden", legend_loc="on data", title="CPIc_C2 high-quality T cells: Leiden clusters", show=False)
    plt.savefig(project_path("results/figures/umap_leiden.png"), dpi=180, bbox_inches="tight")
    plt.close("all")
    logger.info("Transcriptome complete: cells=%d clusters=%d", adata.n_obs, adata.obs["leiden"].nunique())


if __name__ == "__main__":
    main()

