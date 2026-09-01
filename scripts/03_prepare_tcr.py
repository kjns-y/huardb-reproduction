from __future__ import annotations

import argparse

import pandas as pd

from utils import (
    as_bool,
    clean_optional_text,
    load_config,
    normalize_barcode,
    project_path,
    setup_logging,
    validate_tcr_columns,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Select and standardize CPIc_C2 TCR annotations.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    logger = setup_logging("prepare_tcr")

    raw_path = project_path("data/raw") / config["data"]["tcr_file"]
    frame = pd.read_csv(raw_path, compression="gzip", low_memory=False)
    frame = frame.loc[:, ~frame.columns.str.match(r"^Unnamed")].copy()
    validate_tcr_columns(frame)
    suffix = config["barcode"]["tcr_sample_suffix"]
    sample = frame.loc[frame["barcode"].astype(str).str.endswith(suffix)].copy()
    if sample.empty:
        raise ValueError(f"No TCR records found with sample suffix {suffix!r}")

    result = pd.DataFrame({
        "source_barcode": sample["barcode"].astype(str),
        "barcode": sample["barcode"].astype(str).map(
            lambda x: normalize_barcode(x, suffix, config["barcode"]["gex_suffix"])
        ),
        "contig_id": sample["contig_id"].astype(str),
        "chain": clean_optional_text(sample["chain"]),
        "v_gene": clean_optional_text(sample["v_gene"]),
        "d_gene": clean_optional_text(sample["d_gene"]),
        "j_gene": clean_optional_text(sample["j_gene"]),
        "c_gene": clean_optional_text(sample["c_gene"]),
        "cdr3_aa": clean_optional_text(sample["cdr3"]),
        "cdr3_nt": clean_optional_text(sample["cdr3_nt"]).str.upper(),
        "productive": as_bool(sample["productive"]),
        "high_confidence": as_bool(sample["high_confidence"]),
        "is_cell": as_bool(sample["is_cell"]),
        "full_length": as_bool(sample["full_length"]),
        "umis": pd.to_numeric(sample["umis"], errors="coerce").fillna(0).astype(int),
        "reads": pd.to_numeric(sample["reads"], errors="coerce").fillna(0).astype(int),
        "raw_clonotype_id": clean_optional_text(sample["raw_clonotype_id"]),
        "raw_consensus_id": clean_optional_text(sample["raw_consensus_id"]),
    })
    output = project_path("data/processed/tcr/cpic_c2_tcr_contigs.tsv.gz")
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False, compression="gzip")
    sample.to_csv(
        output.parent / "cpic_c2_filtered_contig_annotations.csv.gz",
        index=False,
        compression="gzip",
    )

    summary = pd.DataFrame([{
        "source_rows_all_samples": len(frame),
        "cpic_c2_contigs": len(result),
        "unique_source_barcodes": result["source_barcode"].nunique(),
        "unique_normalized_barcodes": result["barcode"].nunique(),
        "duplicate_contig_ids": int(result["contig_id"].duplicated().sum()),
        "sample_suffix": suffix,
        "target_gex_suffix": config["barcode"]["gex_suffix"],
    }])
    summary.to_csv(project_path("results/qc/tcr_initial_qc.tsv"), sep="\t", index=False)
    logger.info("TCR selected: rows=%d barcodes=%d", len(result), result["barcode"].nunique())


if __name__ == "__main__":
    main()
