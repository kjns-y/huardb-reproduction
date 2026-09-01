# Optional FASTQ → Cell Ranger stage

This stage is deliberately **not run by default**. The processed-data pipeline and final audit must pass first. Cell Ranger is not installed on the server as of 2026-09-01.

## Required SRA runs for CPIc_C2

| Assay | GEO | SRA experiment | Runs | SRA size reported by NCBI |
|---|---|---|---|---:|
| 10x 5′ V1 GEX | GSM4288827 | SRX7647354 | SRR10984983 | 25,301 MB |
| 10x TCR V(D)J | GSM4288865 | SRX7647392 | SRR10985023, SRR10985024 | 2,761 + 2,728 MB |

Source: NCBI SRA RunInfo queried 2026-09-01. The SRA total is about 30.8 GB; paired FASTQ and Cell Ranger intermediates will be substantially larger. Reserve at least 150 GB free, with 250 GB recommended for downloads, FASTQ, references, work directories, and outputs.

## Exact-style software/reference targets

- GEX: Cell Ranger 3.0.2, `cellranger count`, 10x GRCh38 reference v3.0.0.
- TCR: Cell Ranger 3.1.0, `cellranger vdj`, 10x human V(D)J GRCh38 reference v3.1.0.
- The publication reports `cellranger aggr` without depth normalization for the study-wide GEX aggregation. For a single CPIc_C2 reproduction, first validate `count` output against the supplied sample-level processed matrix before considering aggregation.

Modern Cell Ranger versions may change gene annotations, barcode calling, contig filtering, and matrix dimensions. They are a modern conceptual rerun, not an exact software-version reproduction.

## Proposed commands (not executed)

1. Download with SRA Toolkit `prefetch`, then create paired FASTQ with `fasterq-dump --split-files`.
2. Confirm 10x read structure/sample sheet before naming FASTQs for Cell Ranger.
3. Run `cellranger count` for SRR10984983-derived GEX FASTQs.
4. Run `cellranger vdj` for SRR10985023/SRR10985024-derived TCR FASTQs.
5. Compare matrix shape, called barcodes, contig counts, barcode intersection, and strict hcT counts against the processed-data audit.

Large downloads or installation of Cell Ranger require explicit approval and are outside the default pipeline.

