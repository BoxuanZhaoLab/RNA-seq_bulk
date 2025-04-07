# RNA-seq_bulk
 **Input**: genome annotatioin(.gtf), RNA transcriptome(.fasta)

**Output**: alignment(.bam)

**Package used: STAR**

### Building genome index
#### (Required) Download genome(.fa) and annotation(.gtf) from [ncbi genome dataset](https://www.ncbi.nlm.nih.gov/datasets/genome/)

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
