# RNA-seq_bulk

## Run seperately

### A. Preprocessing

#### 1. Generating Genome Index (might already exist, check folder with name '_index')
*Package used: STAR*
**(Required) Download genome(.fa) and annotation(.gtf) from [ncbi genome dataset](https://www.ncbi.nlm.nih.gov/datasets/genome/)**

*Indexing allows the aligner to quickly locate potential matches.*

**Input**: genome annotatioin(.gtf), Reference genome(.fa)

**Output**: index files

- **--runThreadN** `<NumberOfThreads>`  
  Number of CPU threads to use (set to available cores).

- **--runMode** `genomeGenerate`  
  Tells STAR to build genome indices.

- **--genomeDir** `/path/to/genomeDir`  
  Directory (must exist) where indices will be stored. Needs ≥100 GB free for mammalian genomes.

- **--genomeFastaFiles** `/path/to/genome1.fa [/path/to/genome2.fa ...]`  
  One or more reference FASTA files (chromosomes can be split across files).

- **--sjdbGTFfile** `/path/to/annotations.gtf`  
  (Optional) GTF containing exon/junction annotations. Improves splice‐aware mapping.

- **--sjdbOverhang** `<ReadLength − 1>`  
  Length of genomic sequence around each splice junction. Ideally set to ReadLength−1 (e.g. 99 for 100 bp reads).

> **Notes:**  
> - Create (`mkdir`) and clear `genomeDir` before running.  
> - STAR stores indices in binary format—avoid manual editing.  

Example script

```bash
conda activate fastp
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir INDEX_DIR \
--genomeFastaFiles GENOME_DIR/genome.fa \
--sjdbGTFfile GENOME_DIR/genome.gtf \
--sjdbOverhang 99
```

#### 2. rRNA removal
*Package used: ribodetector*

**Input**: RNA sample transctiptome(fastq.gz)

**Output**: no rRNA transcriptome(.fa)

```bash
conda activate ribodetector
./RMrrna.sh /path/to/fastq_dir /path/to/output_dir
```

#### 3. Adapter trimming
*Package used: fastp*

**Input**: no rRNA transcriptome(.fa)

**Output**: trimmed RNA transcriptome(.fastq.gz)

```bash
conda activate fastp
./april_trim.loop.sh
```
> **Notes:**
> Follow prompted instruction to type in the *absolute* directory

### B. Sequence alignment
*Package used: STAR*

**Input**: trimmed RNA transcriptome(.fastq.gz)

**Output**: aligned reads(.bam)

> **Notes:**  
> *INDEX_DIR="~/data/m6A_william/mouseGenome/STAR_mouse_genome_index_ensembl"*
> *INDEX_DIR="~/data/m6A_william/mouseGenome/STAR_mouse_genome_index_ncbi"*

```bash
conda activate fastp
./april_star /path/to/fastq_dir /path/to/output_dir /path/to/index_dir
```

### C. Gene counts
*Package used: HTseq-count*
**Input**: aligned reads(.bam), genome annotatioin(.gtf)

**Output**: count file(.txt)

```bash
conda activate fastp
./april_count /path/to/fastq_dir /path/to/output_dir /path/to/annotation_dir
```

## Run altogether

```bash
./final_version
```
The script will then prompt you for genome species, input directory, and output directory:
```
Enter species (mouse/zebrafish) [mouse]: 
Input directory [$pwd] for trasncriptome raw data: 
Output root directory /path/to/your/output/with/run/date: 
```
> **Notes:** 
> - Create a parent output directory with the command in the following format ```mkdir ~/lab_data/samplename_ddmmyy``` (This should be what you type in for Output root directory)
> - This script will create child directory for each step automatically as "step1_RMrrna", "step2_adapter_trim", "step3_alignment", and "step4_gene_count" To learn about the difference in absolute path and relative path, refrence this link: [Directory(computing)](https://en.wikipedia.org/wiki/Directory_(computing))
> - All directory typed in should be in absolute path format. To learn about the difference in absolute path and relative path, refrence this link: [Path(computing)](https://en.wikipedia.org/wiki/Path_(computing))



