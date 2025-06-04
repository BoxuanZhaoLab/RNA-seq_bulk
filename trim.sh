#!/usr/bin/env bash
set -euo pipefail

# Usage: $0 /path/to/input_dir /path/to/output_root
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 /path/to/input_dir /path/to/output_root" >&2
  exit 1
fi

input_dir="$1"
out_root="$2"
qc_dir="$out_root/qcFiles"

if [[ ! -d "$input_dir" ]]; then
  echo "Error: Input directory '$input_dir' does not exist!" >&2
  exit 1
fi

mkdir -p "$out_root" "$qc_dir"
shopt -s nullglob

for fa in "$input_dir"/*.fa; do
  sample=$(basename "$fa" .fa)

  fastp \
    -i "$fa" \
    -o "$out_root/${sample}.trimmed.fastq.gz" \
    --thread 16 \
    -g \
    -j "$qc_dir/${sample}.json" \
    -h "$qc_dir/${sample}.html" \
    --failed_out "$qc_dir/${sample}.failed.fastq.gz"
done
