import pandas as pd

from scripts.utils import strict_paired_barcodes


def test_strict_chain_filter_keeps_only_one_tra_one_trb() -> None:
    rows = []
    definitions = {
        "cell1": ["TRA", "TRB"],
        "cell2": ["TRA", "TRA", "TRB"],
        "cell3": ["TRA"],
    }
    for barcode, chains in definitions.items():
        for index, chain in enumerate(chains):
            rows.append({
                "barcode": barcode,
                "chain": chain,
                "productive": True,
                "high_confidence": True,
                "is_cell": True,
                "cdr3_nt": f"NT{index}",
            })
    result = strict_paired_barcodes(pd.DataFrame(rows))
    assert result.tolist() == ["cell1"]

