# EpitopePipeline
Pipeline for predicting neoepitopes

## Pipeline description
* The Epitope pipeline is consisted of three parts: peptide generation, hla typing, and neoepitope prediction
### Part 1: Peptide generation

#### 1A. File prep
##### Convert bam to pileup
* samtools 1.3.1

```
samtools mpileup -f fasta bam > pileup
```

##### Index bam
* samtools 1.3.1

```
samtools index bam bai
```

#### 1A. Variant calling
#### 1B. Variant annotation
### Part 2: HLA typing
### Part 3: Neoepitope prediction
