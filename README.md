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
#### Configure directory
- Configure the directory where inputs needed for MHC and outputs from MHC will be stored.
- Use the Python script `configuration.py`:
```
 python configuration.py -h
Usage: python configuration.py <options>

Options:
  -h, --help            show this help message and exit
  --directory=DIRECTORY
  --patientID=PATIENTID
  --peptides_directory=PEPTIDES_DIRECTORY
```

```
python configuration.py --directory {directory_to_create} --peptides_directory {peptides_directory} --patientID {patientID}
```
- The Python script `configuration.py` does the following:
  + create a directory where the inputs for MHC and outputs from MHC will be stored. The name of this directory is specified by the argument `directory`.
  + under the new directory `directory`, create 3 directories: peptides, hla, and IEDB_out. 
  + copy the peptide files to the peptides directory

- Under the directory `hla`, each text file is an HLA type. For example:
```
cat hla_a_03_02_01.txt
HLA-A   hla_a_03_02_01
```

- A note about naming convention, since we are interested in HLA-A, HLA-B, and HLA-C, and there are 2 alleles for each, I named them: hla_a_1, hla_a_2, hla_b_1, hla_b_2, hla_c_1, and hla_c_2. 

#### Configure IEDB tool
1. Download the IEDB tool from http://tools.iedb.org/mhci/download/
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
 
2. Download the file `predict_binding.py` from this Github repo. You will need to replace the `predict_binding.py` script that came with IEDB_MHC with this file.  
- **Important: the tool will not work if you use the script `predict_binding.py` that comes with the installation.**

#### Run MHC
1. Use the script `run_mhc.py`. This Python script has to be run in Python 2. We can use a conda directive to invoke Python2 in snakemake.

```
python run_mhc.py -h
usage: run_mhc.py [-h] --hla HLA --patientID PATIENTID --hla_type HLA_TYPE

Run mhc

optional arguments:
  -h, --help            show this help message and exit
  --hla HLA             REQUIRED. Input the path to the hla file.
  --patientID PATIENTID
                        REQUIRED. Input the patient ID. PatientID is also the
                        name of the output directory
  --hla_type HLA_TYPE   REQUIRED.
```

2. Example:
```
python run_mhc.py --hla TCGA-BC-A10Y_WES_XY_withoutYpar/hla/hla_a_03_02_01.txt --patientID TCGA-BC-A10Y_WES_XY_withoutYpar --hla_type hla_a_1
```

#### For each transcript, merge and select peptide with the lowest IEDB score**
1. Use the script `find_potential_neoepitope.py`. This Python script is run in Python 3.
```
python find_potential_neoepitope.py -h
usage: find_potential_neoepitope.py [-h] --hla HLA --patientID PATIENTID
                                    --hla_type HLA_TYPE

Find potential neoepitope

optional arguments:
  -h, --help            show this help message and exit
  --hla HLA             REQUIRED. Input the path to the hla file.
  --patientID PATIENTID
                        REQUIRED. Input the patient ID. PatientID is also the
                        name of the output directory
  --hla_type HLA_TYPE   REQUIRED.
```

2. Example:
```
python find_potential_neoepitope.py --hla TCGA-BC-A10Y_WES_XY_withoutYpar/hla/hla_a_03_02_01.txt --patientID TCGA-BC-A10Y_WES_XY_withoutYpar --hla_type hla_a_1
```

#### Remove transcript duplicates, and only keep genes where IEDB score <= 500
1. Use the script `filter_neoepitope.py`. This Python script is run in Python 3. 
```
python filter_neoepitope.py -h
Usage: python filter_neoepitope.py <options>

Options:
  -h, --help            show this help message and exit
  --patientID=PATIENTID
  --hla_type=HLA_TYPE
```

2. Example:
```
python filter_neoepitope.py --patientID TCGA-BC-A10Y_WES_XY_withoutYpar --hla_type hla_a_1
```

