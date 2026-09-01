from __future__ import annotations

import argparse

import pandas as pd
import scanpy as sc

from utils import (
    define_huardb_clonotypes,
    load_config,
    project_path,
    setup_logging,
    strict_paired_barcodes,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build strict paired hcT cells and manual huARdb-style clonotypes.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    logger = setup_logging("define_clonotypes")

    tcr = pd.read_csv(project_path("data/processed/tcr/cpic_c2_tcr_contigs.tsv.gz"), sep="\t")
    adata = sc.read_h5ad(project_path("results/objects/cpic_c2_qc.h5ad"), backed="r")
    gex_barcodes = pd.Index(adata.obs_names.astype(str))
    adata.file.close()

    eligible = tcr.loc[
        tcr["productive"] & tcr["high_confidence"] & tcr["is_cell"] & tcr["chain"].isin(["TRA", "TRB"])
    ].copy()
    eligible = eligible.loc[eligible["barcode"].isin(gex_barcodes)].copy()
    strict_barcodes = strict_paired_barcodes(eligible)
    strict_contigs = eligible.loc[eligible["barcode"].isin(strict_barcodes)].copy()

    fields = ["v_gene", "j_gene", "cdr3_aa", "cdr3_nt"]
    tra = strict_contigs.loc[strict_contigs["chain"].eq("TRA"), ["barcode"] + fields].rename(
        columns={field: f"TRA_{field.replace('v_gene', 'v').replace('j_gene', 'j')}" for field in fields}
    )
    trb = strict_contigs.loc[strict_contigs["chain"].eq("TRB"), ["barcode"] + fields].rename(
        columns={field: f"TRB_{field.replace('v_gene', 'v').replace('j_gene', 'j')}" for field in fields}
    )
    hct = tra.merge(trb, on="barcode", how="inner", validate="one_to_one").sort_values("barcode").reset_index(drop=True)
    clones, cell_map = define_huardb_clonotypes(hct)
    expanded_min = int(config["clonotype"]["expanded_min_size"])
    clones["is_expanded"] = clones["clone_size"].ge(expanded_min)
    cell_map["is_expanded"] = cell_map["clone_size"].ge(expanded_min)

    hct_with_clone = hct.merge(
        cell_map[["barcode", "clonotype_id", "clone_size", "clone_rank", "is_expanded"]],
        on="barcode",
        validate="one_to_one",
    )
    hct_with_clone.to_csv(project_path("results/tables/hct_cells.tsv"), sep="\t", index=False)
    clones.to_csv(project_path("results/tables/clonotypes.tsv"), sep="\t", index=False)
    cell_map.to_csv(project_path("results/tables/cell_to_clonotype.tsv"), sep="\t", index=False)

    tcr_barcodes = pd.Index(tcr.loc[tcr["is_cell"], "barcode"].dropna().unique())
    tcr_associated = int(gex_barcodes.intersection(tcr_barcodes).size)
    pairing = pd.DataFrame([{
        "total_gex_cells": len(gex_barcodes),
        "tcr_associated_cells": tcr_associated,
        "strict_hct_cells": len(hct),
        "strict_hct_percentage_of_tcr_associated": 100 * len(hct) / tcr_associated if tcr_associated else 0.0,
        "unique_clonotypes": len(clones),
        "expanded_clonotypes": int(clones["is_expanded"].sum()),
    }])
    pairing.to_csv(project_path("results/qc/paired_tcr_summary.tsv"), sep="\t", index=False)
    logger.info("Strict hcT=%d unique_clonotypes=%d", len(hct), len(clones))


if __name__ == "__main__":
    main()
