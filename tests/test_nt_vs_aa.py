import pandas as pd

from scripts.utils import define_huardb_clonotypes


def test_same_aa_different_nt_are_distinct_clonotypes() -> None:
    frame = pd.DataFrame([
        {"barcode": "cell1", "TRA_v": "TRAV1", "TRA_j": "TRAJ1", "TRA_cdr3_aa": "CASS", "TRA_cdr3_nt": "AAA", "TRB_v": "TRBV1", "TRB_j": "TRBJ1", "TRB_cdr3_aa": "CATS", "TRB_cdr3_nt": "CCC"},
        {"barcode": "cell2", "TRA_v": "TRAV1", "TRA_j": "TRAJ1", "TRA_cdr3_aa": "CASS", "TRA_cdr3_nt": "AAG", "TRB_v": "TRBV1", "TRB_j": "TRBJ1", "TRB_cdr3_aa": "CATS", "TRB_cdr3_nt": "CCC"},
    ])
    _, cells = define_huardb_clonotypes(frame)
    assert cells["clonotype_id"].nunique() == 2

