from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from utils import audit_clone_sizes, project_path, setup_logging


def main() -> None:
    parser = argparse.ArgumentParser(description="Run final huARdb reproduction data constraints audit.")
    parser.parse_args()
    logger = setup_logging("final_audit")
    adata = sc.read_h5ad(project_path("results/objects/cpic_c2_merged.h5ad"), backed="r")
    tcr = pd.read_csv(project_path("data/processed/tcr/cpic_c2_tcr_contigs.tsv.gz"), sep="\t")
    hct = pd.read_csv(project_path("results/tables/hct_cells.tsv"), sep="\t")
    clones = pd.read_csv(project_path("results/tables/clonotypes.tsv"), sep="\t")
    cell_map = pd.read_csv(project_path("results/tables/cell_to_clonotype.tsv"), sep="\t")

    eligible = tcr.loc[tcr["productive"] & tcr["high_confidence"] & tcr["is_cell"] & tcr["chain"].isin(["TRA", "TRB"])]
    counts = eligible.loc[eligible["barcode"].isin(hct["barcode"])].groupby(["barcode", "chain"]).size().unstack(fill_value=0)
    obs_names = pd.Index(adata.obs_names.astype(str))
    obs_map = adata.obs.loc[adata.obs["is_hct"], "clonotype_id"].astype(str).sort_index()
    expected_map = cell_map.set_index("barcode")["clonotype_id"].astype(str).sort_index()

    checks: list[tuple[str, bool, str]] = [
        ("gex_barcode_unique", not obs_names.duplicated().any(), f"n={len(obs_names)}"),
        ("hct_barcode_unique", not hct["barcode"].duplicated().any(), f"n={len(hct)}"),
        ("hct_clonotype_non_null", not hct["clonotype_id"].isna().any(), f"missing={hct['clonotype_id'].isna().sum()}"),
        ("every_hct_exactly_one_TRA", bool((counts.get("TRA", 0) == 1).all()), f"barcodes={len(counts)}"),
        ("every_hct_exactly_one_TRB", bool((counts.get("TRB", 0) == 1).all()), f"barcodes={len(counts)}"),
        ("no_missing_TRA_cdr3_nt", not hct["TRA_cdr3_nt"].isna().any(), f"missing={hct['TRA_cdr3_nt'].isna().sum()}"),
        ("no_missing_TRB_cdr3_nt", not hct["TRB_cdr3_nt"].isna().any(), f"missing={hct['TRB_cdr3_nt'].isna().sum()}"),
        ("clone_rank_consecutive", clones["clone_rank"].tolist() == list(range(1, len(clones) + 1)), f"n={len(clones)}"),
        ("clone_size_nonincreasing", clones["clone_size"].is_monotonic_decreasing, "sorted descending"),
        ("adata_obs_matches_cell_map", obs_map.equals(expected_map), f"adata_hct={len(obs_map)} table_hct={len(expected_map)}"),
    ]
    checks.extend(audit_clone_sizes(clones, cell_map))
    required_files = [
        "results/qc/gex_initial_qc.tsv",
        "results/qc/barcode_merge_audit.tsv",
        "results/tables/clonotypes.tsv",
        "results/tables/clone_by_cluster.tsv",
        "results/objects/cpic_c2_merged.h5ad",
        "results/qc/figure_inventory.tsv",
    ]
    for relative in required_files:
        path = project_path(relative)
        checks.append((f"file_{Path(relative).name}", path.exists() and path.stat().st_size > 0, relative))
    adata.file.close()

    output = project_path("results/qc/final_audit.txt")
    with output.open("w") as handle:
        for name, passed, detail in checks:
            handle.write(f"{'PASS' if passed else 'FAIL'}\t{name}\t{detail}\n")
        handle.write(f"SUMMARY\tpassed={sum(x[1] for x in checks)}\ttotal={len(checks)}\n")
    failures = [name for name, passed, _ in checks if not passed]
    if failures:
        raise SystemExit(f"Final audit failed: {failures}")
    logger.info("Final audit PASS: %d constraints", len(checks))


if __name__ == "__main__":
    main()

