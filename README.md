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
**CAUTION:** Make sure that the naming convention for chromosome is consistent between the bam files and the fasta index file. For instant, if the convention of the bam file is "1", "2", "3", and the convention of the fasta index file is "chr1", "chr2", "chr3", this step will not work. 
To check what the naming convention of the bam file is, you can do:
```
samtools view -H {bam}
```
If your bam file labels chr1 as 1, an easy fix is to modify the fasta index file using `sed`:

```
sed -i 's/chr//g' ucsc.hg19.fasta.fai
```

##### Index bam
* samtools 1.3.1

```
samtools index {bam} {bai}
```

#### 1B. Variant calling
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

#### 1C. Variant annotation
* Variant effect predictor version 86

#### 1D. Generate peptide
* pvacseq version 3.0.5

### Part 2: HLA typing
* Use polysolver
* See https://github.com/tanyaphung/run_hla_polysolver for how to run the software

*Temp notes for running Lynch syndrome samples*

* Working directory is `/home/tphung3/scratch/polysolver`. 

```
PERL5LIB=/home/tphung3/softwares/miniconda3/envs/polysolver/lib/site_perl/5.26.2/ scripts/shell_call_hla_type ~/scratch/LynchSyndrome/bams/MIL582_A2_NL.bam Unknown 1 hg19 STDFQ 0 MIL582_A2_NL
```

### Part 3: Neoepitope prediction

1. Create a directory where the direcoty name is the patient's name. In the example here, I'm using `A7-A26G`.
```
mkdir A7-A26G
```

- We will create 3 sub-directories:

```
mkdir peptides hla IEDB_out
```

2. Run mhc
**PART A: Configure IEDB tool**
* Download the IEDB tool from http://tools.iedb.org/mhci/download/
 - untar the folder:
 ```
 tar -xzvf IEDB_MHC_I-2.19.1.tar.gz
 ```
 - configure the tool and exit out of the directory 
 ```
 cd mhc_i/
 ./configure
 cd ..
 ```
 
 - You will get a message like this after configuring the tool:
 ```
 All prerequisites found!
Copying the standalone-specific netMHCcons template into place
IEDB MHC class I binding prediction tools successfully installed!
Use the command 'python src/predict_binding.py' to get started
 ```
 - Download the file `predict_binding.py` (located in `mhc_i/src/` in this Github repo). You will need to replace the `predict_binding.py` script that came with IEDB_MHC with this file.  

**PART B: run mhc**

- Use the script `run_mhc.py`

```
python run_mhc.py --hla A7-A26G/hla/TCGA-A7-A26G-10.hla.txt --patientID A7-A26G --output_dir A7-A26G/
```

**PART C: For each transcript, select peptide with the lowest IEDB score**

- Use the script `find_potential_neoepitope.py`

```
python find_potential_neoepitope.py --dir A7-A26G --hla A7-A26G/hla/TCGA-A7-A26G-10.hla.txt --patientID A7-A26G
```



