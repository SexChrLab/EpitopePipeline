# Setting up filesnames here:
# Format: normal file first, then tumor
file_ids = ["81056a42-991f-4273-81de-f874806df023", "7892a5c4-f71e-4bf8-b309-251dcf260658"]
file_names = ["e59dcbdcafe075313b13f2c45c1b6df7", "3886522fcb4f50f10124f140ae535107"]

SAMPLES = ["TCGA-A7-A26G-10A-01D-A167-09", "TCGA-A7-A26G-01A-21D-A167-09"]
normal_tumor = ["TCGA-A7-A26G"]
tumor_bam = ["TCGA-A7-A26G-01A-21D-A167-09"]

# path to where bams are
bams_path = "/data/storage/DATASETS/INVESTOR_DATA/BUETOW/BRCA"

# Path to packages not in conda
VarScan_path = "external_scripts/VarScan.v2.3.9.jar"
bamreadcount_path = "external_scripts/bam-readcount"
perl_filter_path = "external_scripts/fpfilter-2.pl"
vep_path = "external_scripts/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl"
plugin_path = "/home/tphung3/.vep/Plugins"
perl5lib_path = "/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/"

rule all:
#     input:
#         expand("bams/{sample}.bam", sample=SAMPLES)
    input:
        expand("peptides/{filename}_VarScan_vep.15.peptide", filename=normal_tumor),
        expand("peptides/{filename}_VarScan_vep.17.peptide", filename=normal_tumor),
        expand("peptides/{filename}_VarScan_vep.19.peptide", filename=normal_tumor),
        expand("peptides/{filename}_VarScan_vep.21.peptide", filename=normal_tumor)

rule mk_symlink_for_bam:
    input:
        expand("/data/storage/DATASETS/INVESTOR_DATA/BUETOW/BRCA/{file_id}/{file_name}.bam", zip, file_id=file_ids, file_name=file_names)
    output:
        expand("bams/{sample}.bam", sample=SAMPLES)
    run:
        for i in range(0,2):
            bam = input[i]
            bam_symlink = output[i]

            shell("ln -s {bam} {bam_symlink}")

rule bam_pileup:
    input:
        fa = "fasta/ucsc.hg19.fasta",
        bam = "bams/{sample}.bam"
    output:
        pileup = "pileups/{sample}.pileup"
    threads: 4
    shell:
        """
        samtools mpileup -f {input.fa} {input.bam} > {output.pileup}
        """

rule index_bam:
    input:
        bam = "bams/{sample}.bam"
    output:
        bai = "bams/{sample}.bam.bai"
    shell:
        """
        samtools index {input.bam} {output.bai}
        """

rule run_VarScan:
    input:
        pileup = expand("pileups/{sample}.pileup", sample=SAMPLES)
    output:
        snp = "intermediate_files/{filename}.VarScan.snp",
        indel = "intermediate_files/{filename}.VarScan.indel"
    params:
        VarScan = VarScan_path,
        basename = "intermediate_files/{filename}.VarScan"
    threads: 4
    shell:
        "java -jar {params.VarScan} somatic {input.pileup} {params.basename} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05"

rule isolate_calls_by_type_and_confidence:
    input:
        VarScan_snp = "intermediate_files/{filename}.VarScan.snp"
    output:
        VarScan_snp = "intermediate_files/{filename}.VarScan.snp.Somatic.hc"
    params:
        VarScan = VarScan_path
    shell:
        """
        java -jar {params.VarScan} processSomatic {input.VarScan_snp}
        """

rule somatic_filter:
    input:
        snp_somatic_hc = "intermediate_files/{filename}.VarScan.snp.Somatic.hc",
        indel = "intermediate_files/{filename}.VarScan.indel"
    output:
        snp_somatic_hc_filter = "intermediate_files/{filename}.VarScan.snp.Somatic.hc.filter"
    params:
        VarScan = VarScan_path
    shell:
        """
        java -jar {params.VarScan} somaticFilter {input.snp_somatic_hc} -indel-file {input.indel} -output-file {output.snp_somatic_hc_filter}
        """

rule convert_to_bed_fmt:
    input:
        snp_somatic_hc_filter = "intermediate_files/{filename}.VarScan.snp.Somatic.hc.filter"
    output:
        snp_somatic_hc_filter_bed = "intermediate_files/{filename}.VarScan.snp.Somatic.hc.filter.bed"
    shell:
        """
        awk -F "\t" '{{print $1 "\t" $2 "\t" $2 }}' {input.snp_somatic_hc_filter} | tail -n+2 > {output.snp_somatic_hc_filter_bed}
        """

rule readcount:
    input:
        fa = "fasta/ucsc.hg19.fasta",
        snp_somatic_hc_filter_bed = "intermediate_files/{filename}.VarScan.snp.Somatic.hc.filter.bed",
        tumor_bam = expand("bams/{tumor}.bam", tumor=tumor_bam),
        bai = expand("bams/{sample}.bam.bai", sample=SAMPLES)
    output:
        readcounts = "intermediate_files/{filename}.readcounts"
    params:
        bamreadcount = bamreadcount_path
    shell:
        """
        {params.bamreadcount} -q 1 -b 20 -f {input.fa} -l {input.snp_somatic_hc_filter_bed} {input.tumor_bam} > {output.readcounts}
        """

rule perl_filter:
    input:
        snp_somatic_hc_filter = "intermediate_files/{filename}.VarScan.snp.Somatic.hc.filter",
        readcounts = "intermediate_files/{filename}.readcounts"
    output:
        out = "intermediate_files/{filename}_VarScan_variants_filter.pass"
    params:
        basename = "intermediate_files/{filename}_VarScan_variants_filter",
        perlfilter = perl_filter_path
    shell:
        """
        perl {params.perlfilter} {input.snp_somatic_hc_filter} {input.readcounts} --output-basename {params.basename}
        """

rule convert_to_vcf:
    input:
        "intermediate_files/{filename}_VarScan_variants_filter.pass"
    output:
        "intermediate_files/{filename}_VarScan_variants_filter.pass.vcf"
    shell:
        """
        python external_scripts/VarScan2_format_converter.py {input} > {output}
        """

rule rm_header_from_vcf:
    input:
        "intermediate_files/{filename}_VarScan_variants_filter.pass.vcf"
    output:
        "intermediate_files/{filename}_VarScan_variants_filter.pass.noheader.vcf"
    shell:
        """
        egrep -v "^#" {input} > {output}
        """

rule run_vep:
    input:
        "intermediate_files/{filename}_VarScan_variants_filter.pass.noheader.vcf"
    output:
        "intermediate_files/{filename}_VarScan_vep.vcf"
    params:
        vep = vep_path,
        plugins = plugin_path,
        perl5lib = perl5lib_path
    shell:
        """
        PERL5LIB={params.perl5lib} perl {params.vep} -i {input} --format vcf --cache --offline --vcf -o {output}  --force_overwrite --plugin Wildtype --dir_plugins {params.plugins} --symbol --terms SO --plugin Downstream
        """

rule pvacseq_convert_vcf:
    input:
        "intermediate_files/{filename}_VarScan_vep.vcf"
    output:
        "intermediate_files/{filename}_VarScan_vep.tsv"
    shell:
        """
        pvacseq convert_vcf {input} {output}
        """

rule generate_fasta:
    input:
        "intermediate_files/{filename}_VarScan_vep.tsv"
    output:
        len_15 = "peptides/{filename}_VarScan_vep.15.peptide",
        len_17 = "peptides/{filename}_VarScan_vep.17.peptide",
        len_19 = "peptides/{filename}_VarScan_vep.19.peptide",
        len_21 = "peptides/{filename}_VarScan_vep.21.peptide",
    shell:
        """
        pvacseq generate_fasta {input} 15 {output.len_15};
        pvacseq generate_fasta {input} 17 {output.len_17};
        pvacseq generate_fasta {input} 19 {output.len_19};
        pvacseq generate_fasta {input} 21 {output.len_21};
        """
