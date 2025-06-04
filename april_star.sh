#!/usr/bin/env bash

DATADIR="$1"
OUTDIR="$(realpath -m "${2:-star_output}")"
INDEX_DIR="$3"
# Check input
if [ -z "$DATADIR" ] || [ ! -d "$DATADIR" ]; then
    echo "Error: Missing or invalid input directory" >&2
    exit 1
fi


# Process files
for f in "$DATADIR"/*.fa; do
    base=$(basename "$f" .fa)
    STAR --runThreadN 10 \
        --genomeDir "$INDEX_DIR" \
        --readFilesIn "$f" \
        --outFileNamePrefix "$OUTDIR/${base}_" \
        --outSAMtype BAM SortedByCoordinate
        # --readFilesCommand zcat \
done