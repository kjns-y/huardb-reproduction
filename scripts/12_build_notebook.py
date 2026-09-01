from __future__ import annotations

from pathlib import Path

import nbformat as nbf
import pandas as pd
from nbclient import NotebookClient

from utils import project_path, setup_logging


def main() -> None:
    logger = setup_logging("build_notebook")
    summary = pd.read_csv(project_path("results/tables/analysis_summary.tsv"), sep="\t").iloc[0]
    top = pd.read_csv(project_path("results/tables/clonotypes.tsv"), sep="\t").iloc[0]
    notebook = nbf.v4.new_notebook()
    notebook["metadata"]["kernelspec"] = {
        "display_name": "Python 3 (huARdb reproduction)",
        "language": "python",
        "name": "python3",
    }
    notebook["metadata"]["language_info"] = {"name": "python", "version": "3.11"}
    notebook["cells"] = [
        nbf.v4.new_markdown_cell(
            "# huARdb CPIc_C2 processed-data reproduction\n\n"
            "## tl;dr\n\n"
            f"After GEX QC, **{int(summary['n_cells']):,}** cells remain; "
            f"**{int(summary['strict_hct_cells']):,}** are strict 1TRA+1TRB hcT cells. "
            f"They form **{int(summary['n_clonotypes']):,}** exact paired CDR3-nt clonotypes, "
            f"including **{int(summary['n_expanded_clones']):,}** expanded clones. "
            f"The largest clone is `{top['clonotype_id']}` with **{int(top['clone_size'])}** cells."
        ),
        nbf.v4.new_markdown_cell(
            "## Context & Methods\n\n"
            "This companion notebook reads pipeline outputs rather than recomputing the full Scanpy workflow. "
            "The executable source of truth is `run_pipeline.sh` and `scripts/`.\n\n"
            "### Key Assumptions\n\n"
            "- TCR barcodes ending in `-C2` are explicitly mapped to the matching GEX `-1` suffix.\n"
            "- Strict hcT means productive, high-confidence, exactly one TRA and one TRB, both with CDR3 nucleotide sequence.\n"
            "- Expanded means clone size ≥2.\n"
            "- Clone-level DEG is exploratory and is not donor-aware inference."
        ),
        nbf.v4.new_markdown_cell("## Data\n\n### 1. Load audited result tables"),
        nbf.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "from IPython.display import display, Image\n\n"
            "ROOT = Path.cwd()\n"
            "summary = pd.read_csv(ROOT / 'results/tables/analysis_summary.tsv', sep='\\t')\n"
            "clones = pd.read_csv(ROOT / 'results/tables/clonotypes.tsv', sep='\\t')\n"
            "scirpy_validation = pd.read_csv(ROOT / 'results/qc/scirpy_clonotype_validation.tsv', sep='\\t')\n"
            "display(summary)"
        ),
        nbf.v4.new_markdown_cell("## Results\n\n### 2. Inspect the largest nucleotide clonotypes"),
        nbf.v4.new_code_cell(
            "columns = ['clonotype_id', 'clone_size', 'TRA_v', 'TRA_cdr3_aa', 'TRA_cdr3_nt', 'TRB_v', 'TRB_cdr3_aa', 'TRB_cdr3_nt']\n"
            "display(clones.loc[:, columns].head(10))"
        ),
        nbf.v4.new_markdown_cell("### 3. Confirm manual-versus-Scirpy validation"),
        nbf.v4.new_code_cell("display(scirpy_validation)\nassert scirpy_validation['match'].all()"),
        nbf.v4.new_markdown_cell("### 4. Display the top-clone transcriptome projection"),
        nbf.v4.new_code_cell("Image(filename=str(ROOT / 'results/figures/umap_top_clones.png'))"),
        nbf.v4.new_markdown_cell(
            "## Takeaways\n\n"
            "- Barcode integration is nearly complete but not perfect; unmatched GEX-only and TCR-only cells remain visible in the merge audit.\n"
            "- Most clonotypes are singletons, while a smaller expanded set creates a long-tailed clone-size distribution.\n"
            "- Expanded clones occupy non-identical Leiden/marker states, supporting coupled clonotype–transcriptome exploration.\n"
            "- A clonotype is a receptor-sequence grouping, not proof of antigen specificity."
        ),
    ]
    output = project_path("notebooks/optional_exploration.ipynb")
    output.parent.mkdir(parents=True, exist_ok=True)
    client = NotebookClient(notebook, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(project_path("."))}})
    executed = client.execute()
    nbf.write(executed, output)
    logger.info("Executed notebook written: %s", output)


if __name__ == "__main__":
    main()

