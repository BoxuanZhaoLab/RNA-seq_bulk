# RNA-seq_bulk

## A. Preprocessing

### 1. Generating Genome Index (might already exist, check folder with name '_index')
#### *Package used: STAR*
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

```
conda activate fastp
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir INDEX_DIR \
--genomeFastaFiles GENOME_DIR/genome.fa \
--sjdbGTFfile GENOME_DIR/genome.gtf \
--sjdbOverhang 99
```

--runThreadN 6:
Sets the number of threads to 6 for parallel processing.

--runMode genomeGenerate:
Specifies that STAR should generate a genome index.

--genomeDir INDEX_DIR:
Defines the directory where the generated genome index output will be stored. *Replace INDEX_DIR with your directory*

--genomeFastaFiles GENOME_DIR/genome.fa:
Provides the path to the genome FASTA file for indexing. *Replace GENOME_DIR with your directory*

### 2. rRNA removal
#### *Package used: ribodetector*

**Input**: RNA sample transctiptome(fastq.gz)

**Output**: no rRNA transcriptome(.fa)

```
conda activate ribodetector
./RMrrna.sh /path/to/fastq_dir /path/to/output_dir
```

### 3. Adapter trimming
#### *Package used: fastp*

**Input**: no rRNA transcriptome(.fa)

**Output**: trimmed RNA transcriptome(.fastq.gz)

```
conda activate fastp
./april_trim.loop.sh
```
> **Notes:**
> Follow prompted instruction to type in the *absolute* directory

## B. Sequence alignment
### *Package used: STAR*

**Input**: trimmed RNA transcriptome(.fastq.gz)

**Output**: aligned reads(.bam)

> **Notes:**  
> *INDEX_DIR="~/data/m6A_william/mouseGenome/STAR_mouse_genome_index_ensembl"*
> *INDEX_DIR="~/data/m6A_william/mouseGenome/STAR_mouse_genome_index_ncbi"*

```
conda activate fastp
./april_star /path/to/fastq_dir /path/to/output_dir /path/to/index_dir
```

## C. Gene counts
### *Package used: HTseq-count*
**Input**: aligned reads(.bam), genome annotatioin(.gtf)

**Output**: count file(.txt)
```
conda activate fastp
./april_count /path/to/fastq_dir /path/to/output_dir /path/to/annotation_dir
```

