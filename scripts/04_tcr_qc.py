from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from utils import chain_configuration, load_config, project_path, setup_logging


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize CPIc_C2 TCR chain quality.")
    parser.add_argument("--config", default=None)
    args = parser.parse_args()
    load_config(args.config)
    logger = setup_logging("tcr_qc")
    frame = pd.read_csv(project_path("data/processed/tcr/cpic_c2_tcr_contigs.tsv.gz"), sep="\t")

    summary = pd.DataFrame([{
        "all_contigs": len(frame),
        "TRA_contigs": int(frame["chain"].eq("TRA").sum()),
        "TRB_contigs": int(frame["chain"].eq("TRB").sum()),
        "productive_TRA": int((frame["chain"].eq("TRA") & frame["productive"]).sum()),
        "productive_TRB": int((frame["chain"].eq("TRB") & frame["productive"]).sum()),
        "high_confidence_productive_TRA": int((frame["chain"].eq("TRA") & frame["productive"] & frame["high_confidence"] & frame["is_cell"]).sum()),
        "high_confidence_productive_TRB": int((frame["chain"].eq("TRB") & frame["productive"] & frame["high_confidence"] & frame["is_cell"]).sum()),
    }])
    summary.to_csv(project_path("results/qc/tcr_summary_qc.tsv"), sep="\t", index=False)

    eligible = frame.loc[
        frame["productive"] & frame["high_confidence"] & frame["is_cell"] & frame["chain"].isin(["TRA", "TRB"])
    ].copy()
    counts = eligible.groupby(["barcode", "chain"]).size().unstack(fill_value=0)
    for chain in ["TRA", "TRB"]:
        if chain not in counts:
            counts[chain] = 0
    counts = counts[["TRA", "TRB"]].rename(columns={"TRA": "n_TRA", "TRB": "n_TRB"}).reset_index()
    counts.to_csv(project_path("results/tables/tcr_chain_count_per_cell.tsv"), sep="\t", index=False)
    counts["configuration"] = [chain_configuration(a, b) for a, b in zip(counts["n_TRA"], counts["n_TRB"])]
    counts.to_csv(project_path("results/tables/tcr_chain_configuration.tsv"), sep="\t", index=False)

    all_productive = frame.loc[frame["productive"] & frame["is_cell"] & frame["chain"].isin(["TRA", "TRB"]), "barcode"].drop_duplicates().sort_values()
    all_productive.to_frame().to_csv(project_path("results/tables/all_productive_t_cells.tsv"), sep="\t", index=False)

    plot_data = counts["configuration"].value_counts().sort_values(ascending=True)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.45 * len(plot_data))))
    ax.barh(plot_data.index, plot_data.values, color="#2A6FBB", edgecolor="#1C2B39")
    ax.set(xlabel="Number of barcodes", ylabel="Productive high-confidence chain configuration", title="CPIc_C2 TCR chain configurations")
    for index, value in enumerate(plot_data.values):
        ax.text(value, index, f" {value}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(project_path("results/figures/tcr_chain_configuration.png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    logger.info("TCR QC: eligible barcodes=%d configurations=%d", len(counts), len(plot_data))


if __name__ == "__main__":
    main()

