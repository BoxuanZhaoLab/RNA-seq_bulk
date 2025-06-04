#!/usr/bin/env bash
set -euo pipefail
source "$HOME/anaconda3/etc/profile.d/conda.sh"
# -------------------------------------------------
# 1. Default genome locations (edit to taste)
# -------------------------------------------------
DEFAULT_MICE="/home/zhouwan2/data/m6A_william/mouseGenome"
DEFAULT_ZF="/home/zhouwan2/data/m6A_william/zebrafishGenome"

# -------------------------------------------------
# 2. Choose species → set STAR index & genome FASTA
# -------------------------------------------------
read -r -p "Enter species (mouse/zebrafish) [mouse]: " species
species=${species:-mouse}
species=${species,,}   # lowercase

case "$species" in
  mouse)
    index="${DEFAULT_MICE}/STAR_mouse_genome_index_ncbi"
    genome="${DEFAULT_MICE}/GCF_000001635.27/ncbi_GCF_000001635.27"
    ;;
  zebrafish|zebra_fish|drerio)
    index="${DEFAULT_ZF}/STAR_zf_genome_index_ncbi"
    genome="${DEFAULT_ZF}/GCF_000002035.6/ncbi_GCF_000002035.6"
    ;;
  *)
    echo "ERROR: unsupported species '$species' (use mouse or zebrafish)" >&2
    exit 1
    ;;
esac

# -------------------------------------------------
# 3. Input / output directories
# -------------------------------------------------
read -r -e -p "Input directory [$(pwd)] for trasncriptome raw data: " in_dir
in_dir=${in_dir:-$(pwd)}

read -r -e -p "Output root directory /path/to/your/output/with/run/date: " out_dir
out_dir=${out_dir:-$(pwd)}
mkdir -p "$out_dir"

# -------------------------------------------------
# 4. Echo final configuration (helps with debugging)
# -------------------------------------------------
cat <<EOF
---------------- Configuration ----------------
Species            : $species
Genome FASTA root  : $genome
STAR index         : $index
Input directory    : $in_dir
Output directory   : $out_dir
------------------------------------------------
EOF

conda activate ribodetector

mkdir "$out_dir/step1_RMrrna"
"/home/zhouwan2/data/m6A_william/all_in_one/RMrrna.sh" "$in_dir" "$out_dir/step1_RMrrna"

conda activate fastp
mkdir "$out_dir/step2_adapter_trim"
"/home/zhouwan2/data/m6A_william/all_in_one/trim.sh" "$out_dir/step1_RMrrna" "$out_dir/step2_adapter_trim"

mkdir "$out_dir/step3_alignment"
"/home/zhouwan2/data/m6A_william/all_in_one/april_star.sh" "$out_dir/step2_adapter_trim" "$out_dir/step3_alignment" "$index"

mkdir "$out_dir/step4_gene_count"
"/home/zhouwan2/data/m6A_william/all_in_one/april_count.sh" "$out_dir/step3_alignment" "$out_dir/step4_gene_count" "$genome/genomic.gtf"