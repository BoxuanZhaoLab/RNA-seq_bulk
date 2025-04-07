#!/bin/bash

# READ ME: to run the script, first run these commands:
# chmod +x process_sra.sh
# ./process_sra.sh

# Define list of SRR IDs (as filenames)
read -p "Enter your SRA filenames separated by spaces: " -a SRR_LIST

# Check if any filenames were provided
if [ "${#SRR_LIST[@]}" -eq 0 ]; then
    echo "Error: No SRA filenames provided!"
    exit 1
fi

# Directory handling for genome, output, gene counts, and GTF file:
read -p "Enter genome index directory path (press Enter to use current directory): " GENOME_DIR
GENOME_DIR=${GENOME_DIR:-$(pwd)}
if [[ ! -d "$GENOME_DIR" ]]; then
    echo "Error: Directory $GENOME_DIR does not exist!"
    exit 1
fi

read -p "Enter output directory path (press Enter to use current directory): " OUTPUT_DIR
OUTPUT_DIR=${OUTPUT_DIR:-$(pwd)}
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Error: Directory $OUTPUT_DIR does not exist!"
    exit 1
fi

read -p "Enter gene counts directory path (press Enter to use current directory): " FILTERED_DIR
FILTERED_DIR=${FILTERED_DIR:-$(pwd)}
if [[ ! -d "$FILTERED_DIR" ]]; then
    echo "Error: Directory $FILTERED_DIR does not exist!"
    exit 1
fi

read -p "Enter genome annotation (GTF) file path: " GTF_FILE
if [[ ! -f "$GTF_FILE" ]]; then
    echo "Error: GTF file $GTF_FILE not found!"
    exit 1
fi

# Ask if prefetching is needed
read -p "Do you want to prefetch and download the SRA? (yes/no): " prefetch_choice

if [[ "$prefetch_choice" == "yes" ]]; then
    # Input prefetch directory path
    read -p "Enter prefetch directory path (press Enter to use current directory): " SRA_DIR
    SRA_DIR=${SRA_DIR:-$(pwd)}
    if [[ ! -d "$SRA_DIR" ]]; then
        echo "Error: Directory $SRA_DIR does not exist!"
        exit 1
    fi
    cd "$SRA_DIR" || { echo "Error: Cannot change to directory $SRA_DIR"; exit 1; }
else
    # Input SRA directory path where FASTQ files are already located
    read -p "Enter SRA directory path: " SRA_DIR
    SRA_DIR=${SRA_DIR:-$(pwd)}
    if [[ ! -d "$SRA_DIR" ]]; then
        echo "Error: Directory $SRA_DIR does not exist!"
        exit 1
    fi
    cd "$SRA_DIR" || { echo "Error: Cannot change to directory $SRA_DIR"; exit 1; }
fi

# Loop through each SRR ID
for SRR_ID in "${SRR_LIST[@]}"; do
    echo "Processing $SRR_ID..."

    if [[ "$prefetch_choice" == "yes" ]]; then
        # Prefetch and convert to FASTQ using fasterq-dump
        echo "Prefetching and dumping FASTQ for $SRR_ID..."
        prefetch "$SRR_ID"
        fasterq-dump "$SRR_ID" -O "$SRA_DIR"
    else
        echo "Skipping prefetch; using existing FASTQ files for $SRR_ID..."
    fi

    # Align with STAR
    echo "Aligning $SRR_ID with STAR..."
    nohup STAR --runThreadN 8 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$SRA_DIR/${SRR_ID}_1.fastq" "$SRA_DIR/${SRR_ID}_2.fastq" \
        --outFileNamePrefix "$OUTPUT_DIR/$SRR_ID" \
        --outSAMtype BAM SortedByCoordinate > "$OUTPUT_DIR/${SRR_ID}.log" 2>&1 &
    
    # Wait for STAR to complete
    wait

    # Run htseq-count
    echo "Running htseq-count for $SRR_ID..."
    nohup htseq-count -f bam -r pos -s no -i gene_id \
        "$OUTPUT_DIR/${SRR_ID}Aligned.sortedByCoord.out.bam" \
        "$GTF_FILE" \
        > "$FILTERED_DIR/${SRR_ID}.txt" 2>&1 &
    
    # Wait for htseq-count to complete
    wait

    # Clean up: remove FASTQ files, STAR output and log if prefetch was used
    if [[ "$prefetch_choice" == "yes" ]]; then
        echo "Cleaning up $SRR_ID..."
        rm -f "$SRA_DIR/${SRR_ID}_1.fastq" "$SRA_DIR/${SRR_ID}_2.fastq"
        rm -f "$OUTPUT_DIR/${SRR_ID}Aligned.sortedByCoord.out.bam"
        rm -f "$OUTPUT_DIR/${SRR_ID}.log"
    fi

    echo "$SRR_ID processing complete."
done

echo "All SRR IDs processed."
