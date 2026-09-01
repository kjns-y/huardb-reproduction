#!/usr/bin/env bash
set -euo pipefail

project="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw="$project/data/raw"
metadata="$project/data/metadata"
extract="$project/data/processed/gex_source"
mkdir -p "$raw" "$metadata" "$extract"

gex_name="GSM4288827_C2-CD3-genes-barcodes-matrix.tar.gz"
tcr_name="GSE144469_TCR_filtered_contig_annotations_all.csv.gz"
gex_url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4288827&file=${gex_name}&format=file"
tcr_url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144469&file=${tcr_name}&format=file"

download_one() {
  local url="$1"
  local output="$2"
  if [[ ! -s "$output" ]]; then
    curl --fail --location --retry 3 --continue-at - --output "${output}.part" "$url"
    mv "${output}.part" "$output"
  fi
}

download_one "$gex_url" "$raw/$gex_name"
download_one "$tcr_url" "$raw/$tcr_name"

if [[ ! -s "$extract/C2-CD3/matrix.mtx.gz" ]]; then
  tar -xzf "$raw/$gex_name" -C "$extract"
fi

manifest="$metadata/data_manifest.tsv"
printf 'accession\tsample_accession\trole\tfile\turl\tdownload_date\tsize_bytes\tmd5\tsha256\n' > "$manifest"
write_manifest_row() {
  local accession="$1" sample_accession="$2" role="$3" file="$4" url="$5"
  local date size md5 sha256
  date="$(stat -c '%y' "$file" | cut -d' ' -f1)"
  size="$(stat -c '%s' "$file")"
  md5="$(md5sum "$file" | cut -d' ' -f1)"
  sha256="$(sha256sum "$file" | cut -d' ' -f1)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$accession" "$sample_accession" "$role" "$(basename "$file")" "$url" \
    "$date" "$size" "$md5" "$sha256" >> "$manifest"
}

write_manifest_row GSM4288827 GSM4288827 GEX "$raw/$gex_name" "$gex_url"
write_manifest_row GSE144469 GSM4288865 TCR "$raw/$tcr_name" "$tcr_url"

test -s "$extract/C2-CD3/features.tsv.gz"
test -s "$extract/C2-CD3/barcodes.tsv.gz"
test -s "$extract/C2-CD3/matrix.mtx.gz"
echo "Processed-data download and manifest: PASS"

