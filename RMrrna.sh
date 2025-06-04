#!/usr/bin/env bash
set -euo pipefail

# Usage: $0 /path/to/fastq_dir /path/to/output_dir
DATADIR="$1"
OUTDIR="${2:-}"

# 1) Check inputs
if [[ -z "$DATADIR" || -z "$OUTDIR" ]]; then
  echo "Usage: $0 /path/to/fastq_dir /path/to/output_dir" >&2
  exit 1
fi

if [[ ! -d "$DATADIR" ]]; then
  echo "Error: Input directory '$DATADIR' does not exist." >&2
  exit 1
fi

# 2) Prepare output directories
mkdir -p "$OUTDIR" "$OUTDIR/logs"

# 3) Process each FASTQ
shopt -s nullglob
for f in "$DATADIR"/*.fastq.gz; do
  base=$(basename "$f")
  base_noext=${base%.fastq.gz}

  # strip off the nucleotide barcode suffix
  prefix=$(echo "$base_noext" | cut -d'_' -f1-2)

  # build output paths
  out_fa="$OUTDIR/${prefix}.nonrrna.fa"
  logf="$OUTDIR/logs/${prefix}.log.txt"

  echo "Processing: $f"
  echo " → FASTA out: $out_fa"
  echo " → LOG out:   $logf"
  echo "--------------------------------------"

  ribodetector -t 20 \
                   -l 100 \
                   -i "$f" \
                   -e rrna \
                   --chunk_size 4000 \
                   -o "$out_fa" \
                   --log "$logf"
done
