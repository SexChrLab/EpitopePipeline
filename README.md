# EpitopePipeline
Pipeline for predicting neoepitopes

## Pipeline description
* The Epitope pipeline is consisted of three parts: peptide generation, hla typing, and neoepitope prediction
### Part 1: Peptide generation

#### 1A. File prep
##### Convert bam to pileup
* samtools 1.3.1

```
samtools mpileup -f {fasta} {bam} > {pileup}
```

##### Index bam
* samtools 1.3.1

```
samtools index {bam} {bai}
```

#### 1A. Variant calling
##### Run VarScan
* VarScan 2.3.9
* Threshold:
  - minimum covererage: 10
  - minimum varian allele frequency: 0.08
  - somatic p value: 0.05
```
java -jar VarScan.v2.3.9.jar somatic {pielup} {out} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05
```

##### Isolate calls by type and confidence
* VarScan 2.3.9

```
java -jar VarScan.v2.3.9.jar processSomatic {VarScan.snp}
```

##### Somatic filter
* VarScan 2.3.9

```
java -jar VarScan.v2.3.9.jar somaticFilter {VarScan.snp.Somatic.hc} -indelfile {VarScan.indel} -output-file {VarScan.snp.Somatic.hc.filter}
```

##### False positive filter
* binary program `bam-readcount` and perl script `fpfilter-2.pl`

#### 1B. Variant annotation
* Variant effect predictor version 86

#### 1C. Generate peptide
* pvacseq version 3.0.5

### Part 2: HLA typing
* Use polysolver
* See https://github.com/tanyaphung/run_hla_polysolver for how to run the software

### Part 3: Neoepitope prediction
