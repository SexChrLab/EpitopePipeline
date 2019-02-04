# EpitopePipeline
Pipeline for predicting neoepitopes

## Pipeline description
* The Epitope pipeline is consisted of three parts: peptide generation, hla typing, and neoepitope prediction
### Part 1: Peptide generation
#### 1A. Variant calling
##### Convert bam to pileup
* samtools 1.3.1

```
samtools mpileup -f fasta bam > pileup
```

#### 1B. Variant annotation
### Part 2: HLA typing
### Part 3: Neoepitope prediction
