#!/bin/bash

# Set default paths
DEFAULT_SPLICE="/home/zhouwan2/data/m6A_william/mouseGenome/hisatIndex2/mouseSplicesites_121924.txt"
DEFAULT_DATA="/home/zhouwan2/data/m6A_william/trimmedData"
DEFAULT_OUTPUT="/home/zhouwan2/data/m6A_william/mappingResults_010525"
DEFAULT_INPUT_SAMPLES="PS19_input PS19_Tau WT_input"
DEFAULT_IP_SAMPLES="PS19_IP WT_IP"

# Set package paths
hisat2="/home/zhouwan2/hisat2/hisat2"
samtools="/home/zhouwan2/samtools-1.21/samtools"

# Get species input
read -p "Enter index species (mouse): " species
if [ "$species" != "mouse" ]; then
    echo "Error: Only mouse is currently supported"
    exit 1
fi
INDEX="/home/zhouwan2/data/m6A_william/mouseGenome/hisatIndex2/mouseGenome"

# Use defaults if Enter pressed
read -p "Press Enter to use default settings or any key to customize"

read -p "Splice path [$DEFAULT_SPLICE]: " SPLICE
SPLICE="${SPLICE:-$DEFAULT_SPLICE}"

read -p "Data path [$DEFAULT_DATA]: " Data
Data="${Data:-$DEFAULT_DATA}"

read -p "Output path [$DEFAULT_OUTPUT]: " Output
Output="${Output:-$DEFAULT_OUTPUT}"
Output=$(realpath -m "$Output")

read -p "Input samples [$DEFAULT_INPUT_SAMPLES]: " input_samples
input_samples="${input_samples:-$DEFAULT_INPUT_SAMPLES}"

read -p "IP samples [$DEFAULT_IP_SAMPLES]: " ip_samples
ip_samples="${ip_samples:-$DEFAULT_IP_SAMPLES}"

# Process samples
process_sample() {
    local sample=$1
    local type=$2
    $hisat2 -x "$INDEX" \
        --known-splicesite-infile "$SPLICE" \
        -k 1 \
        --no-unal \
        --summary-file "${Output}/alignment_summary/${sample}.align_summary" \
        -p 24 \
        -U "${Data}/${sample}.nonrrna.trimmed.fastq.gz" |
    $samtools view -bS |
    $samtools sort -o "${Output}/${sample}.${type}.bam"
}

# Create summary directory
mkdir -p "${Output}/alignment_summary"

# Process all samples
for sample in $input_samples; do
    process_sample "$sample" "input"
done

for sample in $ip_samples; do
    process_sample "$sample" "m6AIP"
done

echo "All alignments completed!"