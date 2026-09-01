from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import scanpy as sc

from utils import load_config, project_path, setup_logging


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge CPIc_C2 GEX and TCR metadata by audited barcode mapping.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    load_config(args.config)
    logger = setup_logging("merge_gex_tcr")

    adata = sc.read_h5ad(project_path("results/objects/cpic_c2_qc.h5ad"))
    tcr = pd.read_csv(project_path("data/processed/tcr/cpic_c2_tcr_contigs.tsv.gz"), sep="\t")
    cell_map = pd.read_csv(project_path("results/tables/cell_to_clonotype.tsv"), sep="\t")

    gex = pd.Index(adata.obs_names.astype(str))
    tcr_barcodes = pd.Index(tcr.loc[tcr["is_cell"], "barcode"].dropna().astype(str).unique())
    if gex.duplicated().any() or tcr_barcodes.duplicated().any():
        raise ValueError("Duplicate barcode detected before GEX/TCR merge")
    mapping = cell_map.set_index("barcode")
    if mapping.index.duplicated().any():
        raise ValueError("cell_to_clonotype has duplicate barcodes")

    adata.obs["has_tcr"] = gex.isin(tcr_barcodes)
    adata.obs["is_hct"] = gex.isin(mapping.index)
    for column in ["clonotype_id", "clone_size", "clone_rank", "is_expanded"]:
        values = mapping[column].reindex(gex)
        adata.obs[column] = values.to_numpy()
    adata.obs["clone_size"] = pd.to_numeric(adata.obs["clone_size"], errors="coerce").astype("Int64")
    adata.obs["clone_rank"] = pd.to_numeric(adata.obs["clone_rank"], errors="coerce").astype("Int64")
    adata.obs["is_expanded"] = adata.obs["is_expanded"].fillna(False).astype(bool)

    audit = pd.DataFrame([{
        "gex_barcodes": len(gex),
        "tcr_barcodes": len(tcr_barcodes),
        "intersection": len(gex.intersection(tcr_barcodes)),
        "gex_only": len(gex.difference(tcr_barcodes)),
        "tcr_only": len(tcr_barcodes.difference(gex)),
        "strict_hct": int(adata.obs["is_hct"].sum()),
        "barcode_mapping": "TCR suffix -C2 explicitly replaced with GEX suffix -1",
    }])
    audit.to_csv(project_path("results/qc/barcode_merge_audit.tsv"), sep="\t", index=False)
    adata.uns["barcode_merge"] = audit.iloc[0].astype(str).to_dict()
    adata.write_h5ad(project_path("results/objects/cpic_c2_merged_preprocessed.h5ad"), compression="gzip")
    logger.info("Barcode merge: GEX=%d intersection=%d strict_hcT=%d", len(gex), audit.loc[0, "intersection"], audit.loc[0, "strict_hct"])


if __name__ == "__main__":
    main()

