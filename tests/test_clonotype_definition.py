import pandas as pd

from scripts.utils import define_huardb_clonotypes


def _frame() -> pd.DataFrame:
    return pd.DataFrame([
        {"barcode": "cell1", "TRA_v": "TRAV1", "TRA_j": "TRAJ1", "TRA_cdr3_aa": "CA", "TRA_cdr3_nt": "AAA", "TRB_v": "TRBV1", "TRB_j": "TRBJ1", "TRB_cdr3_aa": "CB", "TRB_cdr3_nt": "CCC"},
        {"barcode": "cell2", "TRA_v": "TRAV1", "TRA_j": "TRAJ1", "TRA_cdr3_aa": "CA", "TRA_cdr3_nt": "AAA", "TRB_v": "TRBV1", "TRB_j": "TRBJ1", "TRB_cdr3_aa": "CB", "TRB_cdr3_nt": "CCC"},
        {"barcode": "cell3", "TRA_v": "TRAV1", "TRA_j": "TRAJ1", "TRA_cdr3_aa": "CA", "TRA_cdr3_nt": "AAG", "TRB_v": "TRBV1", "TRB_j": "TRBJ1", "TRB_cdr3_aa": "CB", "TRB_cdr3_nt": "CCC"},
    ])


def test_exact_paired_nt_clonotype() -> None:
    _, cells = define_huardb_clonotypes(_frame())
    ids = cells.set_index("barcode")["clonotype_id"]
    assert ids["cell1"] == ids["cell2"]
    assert ids["cell3"] != ids["cell1"]


def test_stable_largest_clone_first() -> None:
    clones, _ = define_huardb_clonotypes(_frame())
    assert clones.iloc[0]["clone_size"] == 2
    assert clones.iloc[0]["clonotype_id"] == "clonotype_000001"

