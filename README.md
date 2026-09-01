# huARdb core analysis reproduction: CPIc_C2

## Background

huARdb links a single cell's transcriptome phenotype to its paired antigen-receptor clonotype. This repository reproduces that core logic for colon CD3+ T cells from checkpoint inhibitor-associated colitis, without reproducing the huARdb website.

The project is based on Wu et al., *huARdb: human Antigen Receptor database for interactive clonotype-transcriptome analysis at the single-cell level* (Nucleic Acids Research, 2022; DOI: 10.1093/nar/gkab857).

## Dataset

- GEO series: [GSE144469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469), “Molecular Pathways of Colon Inflammation Induced by Cancer Immunotherapy”.
- Sample: CPIc_C2 / `C2-CD3` colon CD3+ cells.
- GEX: [GSM4288827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4288827), `GSM4288827_C2-CD3-genes-barcodes-matrix.tar.gz`.
- TCR: [GSM4288865](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4288865), subset from `GSE144469_TCR_filtered_contig_annotations_all.csv.gz`.

URLs, download dates, sizes, MD5 and SHA256 are in `data/metadata/data_manifest.tsv`. Raw public files are never modified.

## Pipeline

```text
GEX
│
├── QC (genes ≥200; mitochondrial fraction ≤20%)
├── normalization + log1p
├── highly variable genes + PCA
├── neighbors + Leiden
└── UMAP
│
│ explicit barcode mapping: TCR -C2 → GEX -1
↓
TCR
│
├── productive + high-confidence + is_cell
├── TRA/TRB chain inventory
├── strict 1TRA + 1TRB pairing
└── exact (TRA CDR3 nt, TRB CDR3 nt) clonotype
│
↓
GEX + TCR
│
├── clonal expansion
├── clone × Leiden / marker state
├── clone × gene expression
├── clone transcriptome correlation + exploratory DEG
└── TRA V × TRB V pairing
```

## Reproduction level

Completed:

1. **Conceptual reproduction** — DONE.
2. **Processed-data reproduction** — DONE with downloaded GEO processed GEX and V(D)J annotations.
3. **Exact software-version reproduction** — NOT DONE; this project uses a pinned modern Python stack.
4. **FASTQ/Cell Ranger reproduction** — NOT RUN; requirements are documented in `scripts/cellranger/README.md`.

## Results

| Metric | Result |
|---|---:|
| Raw GEX cells | 4,209 |
| QC-passed GEX cells | 4,052 |
| TCR-associated QC cells | 3,707 |
| Strict hcT cells | 2,849 |
| Strict hcT / TCR-associated | 76.85% |
| Unique nucleotide clonotypes | 1,587 |
| Singleton clonotypes | 1,253 |
| Expanded clonotypes (size ≥2) | 334 |
| Largest clonotype | clonotype_000001 |
| Largest clone size | 103 cells |
| Largest clone fraction among strict hcT | 3.62% |
| Leiden clusters | 12 |

Manual pandas and Scirpy 0.19 exact-nucleotide clonotype partitions agree for all 2,849 strict hcT cells. This includes clonotype count, singleton count, expanded count, Top10 sizes, and every barcode's clone membership. See `results/qc/scirpy_clonotype_validation.tsv`.

## Run

```bash
cd /mnt/volume2/YEEE/huardb_reproduction
bash run_pipeline.sh all
```

Individual stages:

```bash
bash run_pipeline.sh download
bash run_pipeline.sh gex
bash run_pipeline.sh tcr
bash run_pipeline.sh clonotype
bash run_pipeline.sh merge
bash run_pipeline.sh transcriptome
bash run_pipeline.sh figures
bash run_pipeline.sh notebook
bash run_pipeline.sh test
bash run_pipeline.sh audit
```

The independent environment is `.conda_env`; direct dependencies are pinned in `envs/environment.yml`, with full resolved package records in `envs/pip-freeze.txt` and `envs/conda-explicit.txt`. Every stage writes a timestamped log under `logs/` and stops on error.

## Key outputs

- `results/objects/cpic_c2_merged.h5ad`: QC-passed, log-normalized GEX with TCR/clone metadata, PCA, Leiden and UMAP.
- `results/objects/cpic_c2_receptor_scirpy.h5ad`: Scirpy AIRR representation for strict hcT cells.
- `results/tables/hct_cells.tsv`: one paired TRA/TRB row per strict hcT cell.
- `results/tables/clonotypes.tsv`: stable nucleotide clonotypes and clone sizes.
- `results/tables/cell_to_clonotype.tsv`: cell-to-clone membership.
- `results/tables/analysis_summary.tsv`: final metrics and diversity extensions.
- `results/figures/`: QC, marker, clonotype, correlation, DEG and V-pairing figures.
- `results/qc/final_audit.txt`: final invariant checks.
- `notebooks/optional_exploration.ipynb`: executed companion notebook reading audited outputs.

## Deviations from the paper

| Paper / source-study workflow | This project |
|---|---|
| Cell Ranger 3.0.2 GEX / 3.1.0 V(D)J | Uses published processed outputs; no FASTQ rerun |
| Cells with DoubletFinder high scores removed | Doublet removal not run; explicitly recorded in `gex_filtering_audit.tsv` |
| SingleR/Monaco cell-subtype prediction | Leiden plus exploratory multi-gene marker scores; no claim of exact subtype reproduction |
| tSNE-centered presentation | UMAP is the primary modern embedding |
| Historical R/software stack | Pinned Python 3.11 + Scanpy 1.10.4 + Scirpy 0.19.0 |

The source GEX matrix has 24 duplicated gene symbols. Original symbols are retained in `adata.var['gene_symbol']`; AnnData index names are made unique for safe access. The TCR all-sample file contains a non-standard `Multi` locus row, which Scirpy warns about during import; the strict TRA/TRB subset is unaffected.

Dual-TRA and other multi-chain cells are retained in processed TCR tables and chain-configuration results. The strict 1TRA+1TRB subset is an engineering definition, not a claim that all excluded cells are artifacts.

## Statistical interpretation

Clone transcriptome correlation uses mean log-normalized expression across highly variable genes for clones with size ≥5 and is exploratory. Top-clone DEG uses single cells and therefore has pseudoreplication risk. Cross-patient inference should use donor-aware pseudobulk or mixed models.

Clonotype does not imply known antigen specificity. Expansion alone does not establish tumor specificity, colitis causality, or any particular antigen.

## Tests and audit

`pytest` covers barcode intersection, explicit suffix mapping, strict chain filtering, paired nucleotide clonotype grouping, stable ID assignment, and the critical “same amino acid / different nucleotide” case.

The final audit fails non-zero unless barcodes are unique, every strict hcT has exactly one TRA and TRB with nucleotide CDR3, clone sizes match membership, clone ranks are stable, and AnnData metadata matches the tabular cell map.

## What I learned from reproducing huARdb

1. **How are scRNA-seq and scTCR-seq linked?** By cell barcode, but only after an explicit sample-aware mapping. In these files GEX uses `-1` while the aggregated TCR table uses `-C2`; blindly stripping suffixes would make provenance ambiguous.
2. **How is paired TRA/TRB defined?** Here it requires a called cell, productive/high-confidence receptor, exactly one TRA and one TRB, and valid CDR3 nucleotide for both chains.
3. **Why did huARdb use 1TRA+1TRB?** It creates an unambiguous one-row receptor and stable clone key for database-scale analysis. The trade-off is exclusion of biologically plausible dual-α cells.
4. **Why use CDR3 nucleotide?** It preserves rearrangement identity. This run's unit test demonstrates that equal CDR3 amino acids encoded by different nucleotides remain different huARdb-style clonotypes.
5. **What does clone size mean?** The number of strict hcT cells sharing the exact paired TRA/TRB nucleotide key in this sample after GEX QC. It is not a patient-level replicate count.
6. **How are expanded clones mapped to T-cell state?** `clonotype_id` is written into `adata.obs`, so the same UMAP/Leiden and marker-expression coordinates can highlight each clone and form clone×cluster/expression tables.
7. **Clonotype versus antigen specificity?** A clonotype is a sequence-defined cell group. Antigen specificity needs external annotation or experimental evidence; it cannot be inferred solely from expansion.
8. **Biggest limitation?** Strict pairing and single-sample cell-level analysis simplify complex receptor biology and cannot support donor-level causal inference. Marker annotation here is exploratory rather than a reference-mapped exact reproduction.
9. **How should huARdb be rebuilt today?** Keep raw AIRR/MuData-compatible receptor representations, preserve dual-chain alternatives, version every barcode transform, use donor-aware statistics, expose exact provenance, and validate manual definitions against a receptor toolkit.
10. **How can this scale to a T-cell atlas?** Parameterize sample-specific suffix maps and QC, retain donor/study/batch columns, process per sample into standardized AnnData/AIRR objects, use stable sequence hashes plus study-scoped IDs, and run cohort-level donor-aware models after per-sample audits pass.

