#!/usr/bin/env bash
set -euo pipefail

project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
env_prefix="$project/.conda_env"
python="$env_prefix/bin/python"
conda_bin="/home/jinsq/anaconda3/bin/conda"
mkdir -p "$project/logs"

ensure_env() {
  if [[ ! -x "$python" ]]; then
    export CONDA_PKGS_DIRS="$project/.conda_pkgs"
    "$conda_bin" env create --prefix "$env_prefix" --file "$project/envs/environment.yml"
  fi
}

run_logged() {
  local label="$1"
  shift
  local stamp
  stamp="$(date +%Y%m%d_%H%M%S)"
  "$@" 2>&1 | tee "$project/logs/${label}_${stamp}.log"
}

run_stage() {
  local stage="$1"
  case "$stage" in
    env)
      ensure_env
      run_logged env "$python" -c 'import sys, scanpy, scirpy; print(sys.version); print("scanpy", scanpy.__version__); print("scirpy", scirpy.__version__)'
      ;;
    download)
      run_logged download bash "$project/scripts/00_download_data.sh"
      ;;
    gex)
      ensure_env
      run_logged gex_prepare "$python" "$project/scripts/01_prepare_gex.py"
      run_logged gex_qc "$python" "$project/scripts/02_gex_qc.py"
      ;;
    tcr)
      ensure_env
      run_logged tcr_prepare "$python" "$project/scripts/03_prepare_tcr.py"
      run_logged tcr_qc "$python" "$project/scripts/04_tcr_qc.py"
      ;;
    clonotype)
      ensure_env
      run_logged clonotype_manual "$python" "$project/scripts/05_define_clonotypes.py"
      if [[ -f "$project/scripts/05b_validate_scirpy.py" ]]; then
        run_logged clonotype_scirpy "$python" "$project/scripts/05b_validate_scirpy.py"
      fi
      ;;
    merge)
      ensure_env
      run_logged merge "$python" "$project/scripts/06_merge_gex_tcr.py"
      ;;
    transcriptome)
      ensure_env
      run_logged transcriptome "$python" "$project/scripts/07_transcriptome_analysis.py"
      ;;
    figures)
      ensure_env
      run_logged clone_analysis "$python" "$project/scripts/08_clonotype_analysis.py"
      run_logged clone_phenotype "$python" "$project/scripts/09_clone_phenotype_analysis.py"
      run_logged figures "$python" "$project/scripts/10_make_figures.py"
      ;;
    notebook)
      ensure_env
      run_logged notebook "$python" "$project/scripts/12_build_notebook.py"
      ;;
    test|tests)
      ensure_env
      run_logged tests "$python" -m pytest -q "$project/tests"
      ;;
    audit)
      ensure_env
      run_logged audit "$python" "$project/scripts/11_final_audit.py"
      ;;
    all)
      run_stage env
      for item in download gex tcr clonotype merge transcriptome figures notebook test audit; do
        run_stage "$item"
      done
      ;;
    *)
      echo "Usage: bash run_pipeline.sh {all|env|download|gex|tcr|clonotype|merge|transcriptome|figures|notebook|test|audit}" >&2
      exit 2
      ;;
  esac
}

run_stage "${1:-all}"
