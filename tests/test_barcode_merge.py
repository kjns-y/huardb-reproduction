from scripts.utils import normalize_barcode


def test_barcode_intersection() -> None:
    gex = {"A", "B", "C"}
    tcr = {"B", "C", "D"}
    assert gex & tcr == {"B", "C"}


def test_explicit_c2_to_gex_suffix_mapping() -> None:
    assert normalize_barcode("AAACCC-C2", "-C2", "-1") == "AAACCC-1"


def test_unexpected_suffix_fails() -> None:
    try:
        normalize_barcode("AAACCC-NC2", "-C2", "-1")
    except ValueError:
        pass
    else:
        raise AssertionError("Unexpected sample suffix must not be silently stripped")

