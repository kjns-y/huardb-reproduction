from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_config(path: str | Path | None = None) -> dict:
    config_path = Path(path) if path else PROJECT_ROOT / "config" / "config.yaml"
    with config_path.open() as handle:
        return yaml.safe_load(handle)


def project_path(relative: str | Path) -> Path:
    return PROJECT_ROOT / relative


def setup_logging(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def as_bool(series: pd.Series) -> pd.Series:
    truthy = {"true", "t", "1", "yes"}
    falsy = {"false", "f", "0", "no"}

    def convert(value: object) -> bool:
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        text = str(value).strip().lower()
        if text in truthy:
            return True
        if text in falsy or text in {"none", "nan", ""}:
            return False
        raise ValueError(f"Cannot parse boolean value: {value!r}")

    return series.map(convert).astype(bool)


def clean_optional_text(series: pd.Series) -> pd.Series:
    result = series.astype("string").str.strip()
    return result.mask(result.str.lower().isin({"", "none", "nan", "na"}), pd.NA)


def normalize_barcode(
    barcode: str,
    source_sample_suffix: str = "-C2",
    target_suffix: str = "-1",
) -> str:
    """Map one known aggregated TCR sample suffix to its GEX library suffix.

    This intentionally does not strip arbitrary suffixes. Unexpected barcodes fail
    loudly so that barcode changes remain auditable.
    """
    if not isinstance(barcode, str) or not barcode:
        raise ValueError("barcode must be a non-empty string")
    if not barcode.endswith(source_sample_suffix):
        raise ValueError(
            f"Barcode {barcode!r} does not end with expected suffix "
            f"{source_sample_suffix!r}"
        )
    core = barcode[: -len(source_sample_suffix)]
    if not re.fullmatch(r"[ACGTN]+", core):
        raise ValueError(f"Unexpected 10x barcode core: {core!r}")
    return f"{core}{target_suffix}"


def validate_tcr_columns(frame: pd.DataFrame) -> None:
    required = {
        "barcode",
        "chain",
        "v_gene",
        "d_gene",
        "j_gene",
        "c_gene",
        "cdr3",
        "cdr3_nt",
        "productive",
        "umis",
        "reads",
        "high_confidence",
        "is_cell",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"TCR file is missing required columns: {missing}")


def chain_configuration(n_tra: int, n_trb: int) -> str:
    if n_tra <= 2 and n_trb <= 2 and not (n_tra == 2 and n_trb == 2):
        return f"{n_tra}TRA_{n_trb}TRB"
    return "multi_chain"


def strict_paired_barcodes(frame: pd.DataFrame) -> pd.Index:
    """Return barcodes with exactly one valid high-confidence TRA and one TRB."""
    required = {"barcode", "chain", "productive", "high_confidence", "is_cell", "cdr3_nt"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Missing strict-pairing columns: {sorted(missing)}")
    eligible = frame.loc[
        frame["productive"]
        & frame["high_confidence"]
        & frame["is_cell"]
        & frame["chain"].isin(["TRA", "TRB"])
    ].copy()
    counts = eligible.groupby(["barcode", "chain"]).size().unstack(fill_value=0)
    for chain in ["TRA", "TRB"]:
        if chain not in counts:
            counts[chain] = 0
    paired = counts.index[(counts["TRA"] == 1) & (counts["TRB"] == 1)]
    candidate = eligible.loc[eligible["barcode"].isin(paired)]
    has_nt = candidate.groupby("barcode")["cdr3_nt"].apply(lambda values: values.notna().all())
    return pd.Index(has_nt.index[has_nt]).sort_values()


def define_huardb_clonotypes(hct: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assign stable clone IDs using exact paired TRA/TRB CDR3 nucleotide keys."""
    required = {"barcode", "TRA_cdr3_nt", "TRB_cdr3_nt"}
    missing = required - set(hct.columns)
    if missing:
        raise ValueError(f"Missing clonotype columns: {sorted(missing)}")
    if hct["barcode"].duplicated().any():
        raise ValueError("strict hcT table contains duplicate barcodes")
    if hct[["TRA_cdr3_nt", "TRB_cdr3_nt"]].isna().any().any():
        raise ValueError("strict hcT table contains missing CDR3 nucleotide sequence")

    working = hct.copy()
    keys = ["TRA_cdr3_nt", "TRB_cdr3_nt"]
    sizes = working.groupby(keys, dropna=False).size().rename("clone_size").reset_index()
    sizes = sizes.sort_values(
        ["clone_size", "TRA_cdr3_nt", "TRB_cdr3_nt"],
        ascending=[False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    sizes["clone_rank"] = np.arange(1, len(sizes) + 1, dtype=int)
    sizes["clonotype_id"] = sizes["clone_rank"].map(lambda x: f"clonotype_{x:06d}")

    cell_map = working.merge(sizes, on=keys, how="left", validate="many_to_one")
    cell_map["is_expanded"] = cell_map["clone_size"].ge(2)

    receptor_columns = [
        "TRA_v",
        "TRA_j",
        "TRA_cdr3_aa",
        "TRA_cdr3_nt",
        "TRB_v",
        "TRB_j",
        "TRB_cdr3_aa",
        "TRB_cdr3_nt",
    ]
    representative_columns = [
        "TRA_v",
        "TRA_j",
        "TRA_cdr3_aa",
        "TRB_v",
        "TRB_j",
        "TRB_cdr3_aa",
    ]
    representatives = (
        cell_map.sort_values("barcode", kind="mergesort")
        .groupby("clonotype_id", sort=False)[representative_columns]
        .first()
        .reset_index()
    )
    clones = sizes.merge(representatives, on="clonotype_id", validate="one_to_one")
    clones["is_expanded"] = clones["clone_size"].ge(2)
    clones = clones[
        ["clonotype_id", "clone_size", "clone_rank", "is_expanded"] + receptor_columns
    ]
    return clones, cell_map


def audit_clone_sizes(clones: pd.DataFrame, cell_map: pd.DataFrame) -> list[tuple[str, bool, str]]:
    observed = cell_map.groupby("clonotype_id").size().sort_index()
    expected = clones.set_index("clonotype_id")["clone_size"].sort_index()
    checks = [
        (
            "clone_size_matches_membership",
            observed.equals(expected.astype(observed.dtype)),
            f"observed_clones={len(observed)} expected_clones={len(expected)}",
        ),
        (
            "sum_clone_size_equals_hct",
            int(clones["clone_size"].sum()) == len(cell_map),
            f"sum={int(clones['clone_size'].sum())} hct={len(cell_map)}",
        ),
    ]
    return checks


def available_gene_symbols(adata: object, genes: Iterable[str]) -> list[str]:
    names = set(getattr(adata, "var_names"))
    return [gene for gene in genes if gene in names]
