from __future__ import annotations

import argparse

import pandas as pd
import scirpy as ir

from utils import load_config, normalize_barcode, project_path, setup_logging


def membership_signature(labels: pd.Series) -> pd.Series:
    groups = labels.groupby(labels).apply(lambda values: tuple(sorted(values.index.astype(str))))
    return labels.map(groups)


def clone_metrics(labels: pd.Series) -> dict[str, int | str]:
    sizes = labels.value_counts()
    return {
        "n_clonotypes": int(len(sizes)),
        "n_singletons": int(sizes.eq(1).sum()),
        "n_expanded": int(sizes.ge(2).sum()),
        "largest_clone_size": int(sizes.max()),
        "top10_sizes": ",".join(map(str, sizes.sort_values(ascending=False).head(10).tolist())),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate manual nucleotide clonotypes with Scirpy 0.19.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    logger = setup_logging("validate_scirpy")

    input_path = project_path("data/processed/tcr/cpic_c2_filtered_contig_annotations.csv.gz")
    receptor = ir.io.read_10x_vdj(input_path, filtered=True)
    suffix = config["barcode"]["tcr_sample_suffix"]
    receptor.obs_names = [
        normalize_barcode(str(barcode), suffix, config["barcode"]["gex_suffix"])
        for barcode in receptor.obs_names
    ]
    manual = pd.read_csv(project_path("results/tables/cell_to_clonotype.tsv"), sep="\t").set_index("barcode")
    receptor = receptor[receptor.obs_names.isin(manual.index)].copy()
    if receptor.n_obs != len(manual):
        missing = manual.index.difference(receptor.obs_names)
        raise ValueError(f"Scirpy receptor object is missing {len(missing)} strict hcT barcodes")

    ir.pp.index_chains(receptor)
    ir.tl.chain_qc(receptor)
    ir.tl.define_clonotypes(
        receptor,
        key_added="scirpy_clone_id",
        receptor_arms="all",
        dual_ir="primary_only",
        within_group=None,
    )
    scirpy_labels = receptor.obs["scirpy_clone_id"].astype("string")
    if scirpy_labels.isna().any():
        raise ValueError(f"Scirpy left {int(scirpy_labels.isna().sum())} strict hcT cells unassigned")
    manual_labels = manual.loc[receptor.obs_names, "clonotype_id"].astype("string")
    manual_labels.index = receptor.obs_names

    manual_signature = membership_signature(manual_labels)
    scirpy_signature = membership_signature(scirpy_labels)
    membership_equal = manual_signature.eq(scirpy_signature)
    comparison = pd.DataFrame({
        "manual_clonotype_id": manual_labels,
        "scirpy_clonotype_id": scirpy_labels,
        "same_membership": membership_equal,
        "chain_pairing": receptor.obs["chain_pairing"].astype(str),
    })
    comparison.index.name = "barcode"
    comparison.to_csv(project_path("results/tables/manual_vs_scirpy_membership.tsv"), sep="\t")

    manual_metrics = clone_metrics(manual_labels)
    scirpy_metrics = clone_metrics(scirpy_labels)
    rows = []
    for metric in manual_metrics:
        rows.append({
            "metric": metric,
            "manual": manual_metrics[metric],
            "scirpy": scirpy_metrics[metric],
            "match": manual_metrics[metric] == scirpy_metrics[metric],
        })
    rows.append({
        "metric": "barcode_membership",
        "manual": len(manual_labels),
        "scirpy": int(membership_equal.sum()),
        "match": bool(membership_equal.all()),
    })
    validation = pd.DataFrame(rows)
    validation.to_csv(project_path("results/qc/scirpy_clonotype_validation.tsv"), sep="\t", index=False)
    receptor.write_h5ad(project_path("results/objects/cpic_c2_receptor_scirpy.h5ad"), compression="gzip")
    if not validation["match"].all():
        logger.warning("Manual and Scirpy clonotypes differ; see validation table. No results were forced to match.")
    else:
        logger.info("Manual and Scirpy exact-NT clonotype partitions match for all %d strict hcT cells", len(manual_labels))


if __name__ == "__main__":
    main()

