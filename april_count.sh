#!/bin/bash

input_dir="$1"
output_dir="$2"
gtf_file="$3"


# Process all BAM files
for bam_file in "$input_dir"/*.bam; do
    sample=$(basename "$bam_file" .bam)
    htseq-count -f bam -r pos -s no -i gene_id "$bam_file" "$gtf_file" > "${output_dir}/${sample}_counts.txt"
done

echo "Gene counting complete! Results in: $output_dir"