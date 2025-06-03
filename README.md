# RNA-seq_bulk

## Preprocessing

### Generating Genome Index 
### *Package used: STAR*
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

--sjdbGTFfile GENOME_DIR/genome.gtf:
Supplies the GTF file containing gene annotations for identifying splice junctions. *Replace GENOME_DIR with your directory*

--sjdbOverhang 99:
Sets the length of the genomic sequence (usually read length minus 1) used to improve splice junction detection. 
