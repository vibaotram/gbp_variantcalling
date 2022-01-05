import os
from Bio import SeqIO
import re

configfile: "config.yaml"

## input
fastq_dir = config["fastq_dir"]
sample,fq, = glob_wildcards(os.path.join(fastq_dir, "{sample}/{fq, .*.(fq|fq.gz|fastq|fastq.gz)}"))
fastq_files, = glob_wildcards(os.path.join(fastq_dir, "{fastq_files, .*/.*.(fq|fq.gz|fastq|fastq.gz)}"))

ref = config["ref"]
chrom = []
fasta_sequences = SeqIO.parse(open(ref),'fasta')
for fasta in fasta_sequences:
    seqname = fasta.id
    if not re.match("Contig", seqname):
        chrom.append(seqname)

## output
# chrom = {"Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11"}

output_dir = config["output_dir"]

gatk_VF_opt = config["gatk_VariantFiltration"]["params"]
gatk_SV_opt = config["gatk_SelectVariants"]["params"]
vcftools_opt = config["vcftools"]["params"]

gatk_VF_suf = "" # config["gatk_VariantFiltration"]["suffix"]
gatk_SV_suf = "" # config["gatk_SelectVariants"]["suffix"]
vcftools_suf = "" # config["vcftools"]["suffix"]

filtered = ""
if gatk_VF_opt:
    gatk_VF_suf = config["gatk_VariantFiltration"]["suffix"]
    filtered += gatk_VF_suf
if gatk_SV_opt:
    gatk_SV_suf = config["gatk_SelectVariants"]["suffix"]
    filtered += gatk_SV_suf
if vcftools_opt:
    vcftools_suf = config["vcftools"]["suffix"]
    filtered += vcftools_suf


if gatk_VF_opt or gatk_SV_opt or vcftools_opt:
    final_dir = "filtered_vcf_by_chrom"
else:
    final_dir = "vcf_by_chrom"


container: "docker://continuumio/miniconda3:4.4.10"

## quality control before processing
rule fastqc:
    input: expand(os.path.join(output_dir, "fastqc/{fastq_files}_fastqc.html"), fastq_files = fastq_files)

rule fastqc_ind:
    input: os.path.join(fastq_dir, "{fastq_files}")
    output:
        html = os.path.join(output_dir, "fastqc/{fastq_files}_fastqc.html"),
        zip = os.path.join(output_dir, "fastqc/{fastq_files}_fastqc.zip"),
    params:
        opt = config["fastqc"]["params"]
    threads: config["fastqc"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        outdir=$(dirname {output.html})
        mkdir $outdir
        fastqc -o $outdir -f fastq
        """


rule gbp_variantcalling:
    input:
        vcf_gz = expand(os.path.join(output_dir, "{final_dir}/all_final{filtered}.vcf.gz"), final_dir = final_dir, filtered = filtered),
        vcf_tbi = expand(os.path.join(output_dir, "{final_dir}/all_final{filtered}.vcf.gz.tbi"), final_dir = final_dir, filtered = filtered),
        # expand(os.path.join(output_dir, "{final_dir}/{chrom}{filtered}.vcf"), final_dir = final_dir, chrom = chrom, filtered = filtered)

## pre-processing if neccessary

# rule fastp

## variant calling from fastq
rule index_ref:
    input: ref
    output: multiext(ref, ".bwt", ".pac", ".ann", ".amb", ".sa") # bwt. pac, ann, amb, sa
    shell:
        """
        bwa index {ref}
        """


rule bwa_mem:
    input:
        fastq = os.path.join(fastq_dir, "{sample}"),
        # ref_index = rules.index_ref.output
    output: temp(os.path.join(output_dir, "bwa_mem/{sample}.sam"))
    params: config["bwa_mem"]["params"]
    threads: config["bwa_mem"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        fastq=$(find {input.fastq} -regex ".*\(1\|2\).\(fq.gz\|fastq.gz\|fq\|fastq\)" | xargs -n1 | sort | xargs)
        echo "bwa mem -t {threads} {params} {ref} $fastq > {output}"
        bwa mem {ref} $fastq -t {threads} {params} -R $(./read_group.sh $fastq)> {output}
        """

rule picard_SortSam:
    input: rules.bwa_mem.output
    output:
        bam = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}.bam")),
        idx = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}.bam.bai")),
    params: config["picard_SortSam"]["params"]
    conda: "conda.yaml"
    shell:
        """
        picard SortSam -I {input} -O {output.bam} -CREATE_INDEX true {params}
        """

rule samtools_view:
    input: rules.picard_SortSam.output.bam
    output:
        bam = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}_samtoolsView.bam")),
        idx = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}_samtoolsView.bam.bai")),
    params: config["samtools_view"]["params"]
    threads: config["samtools_view"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        samtools view {params} -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        """

rule MarkDuplicates:
    input: rules.samtools_view.output.bam
    output:
        bam = os.path.join(output_dir, "sorted_bam/{sample}/{sample}_MarkDuplicate.bam"),
    params: config["MarkDuplicates"]["params"]
    threads: config["MarkDuplicates"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        ref_name=$(basename -s .fa {ref})
        ref_dict="$(dirname {ref})/$ref_name.dict"
        if [[ ! -f $ref_dict ]]
        then
            picard CreateSequenceDictionary -REFERENCE {ref} -OUTPUT $ref_dict
        fi
        gatk MarkDuplicatesSpark -R {ref} -I {input} -O {output} --spark-runner LOCAL
        """

rule HaplotypeCaller:
    input: rules.MarkDuplicates.output
    output: temp(os.path.join(output_dir, "HaplotypeCaller/{sample}.g.vcf")),
    params: config["HaplotypeCaller"]["params"]
    threads: config["HaplotypeCaller"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        if [[ ! -f {ref}.fai ]]
        then
            samtools faidx {ref} --fai-idx {ref}.fai
        fi
        gatk HaplotypeCaller -R {ref} -I {input} -O {output} --native-pair-hmm-threads {threads} -ERC GVCF {params}
        """

GenomicDBImport_input = ' -V '.join(expand(rules.HaplotypeCaller.output, sample = set(sample)))

rule GenomicDBImport:
    input: expand(rules.HaplotypeCaller.output, sample = set(sample))
    output: temp(directory(os.path.join(output_dir, "gendb/{chrom}")))
    params: config["GenomicDBImport"]["params"]
    threads: config["GenomicDBImport"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        gatk GenomicsDBImport -V {GenomicDBImport_input} --genomicsdb-workspace-path {output} --intervals {wildcards.chrom} {params} --reader-threads {threads}
        """

rule CombineGVCFs:
    input: rules.GenomicDBImport.output
    output: os.path.join(output_dir, "vcf_by_chrom/{chrom}.vcf")
    params: config["CombineGVCFs"]["params"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        ref=$(realpath {ref})
        output=$(realpath {output})
        input=$(basename {input})
        cd $(dirname {input})
        gatk CombineGVCFs -R $ref -V gendb://$input -O $output {params}
        """



## variant filtration

x = [gatk_VF_suf, gatk_SV_suf, vcftools_suf]

rule filter_variants:
    input: rules.CombineGVCFs.output
    output: os.path.join(output_dir, "filtered_vcf_by_chrom/{chrom}{filtered}.vcf".format(chrom = "{chrom}", filtered = filtered))
    params:
        gatk_VF_opt = gatk_VF_opt,
        gatk_SV_opt = gatk_SV_opt,
        vcftools_opt = vcftools_opt,
        gatk_VF_suf = gatk_VF_suf,
        gatk_SV_suf = gatk_SV_suf,
        vcftools_suf = vcftools_suf,
    conda: "conda.yaml"
    shell:
        """
        in_name=$(basename {input} | sed 's/.vcf//')
        dir=$(dirname {output})
        input="{input}"
        gatk_VF_opt="{params.gatk_VF_opt}"
        if [[ ! -z $gatk_VF_opt ]]
        then
            echo -e "### filtering by gatk VariantFiltration: {params.gatk_VF_opt}"
            gatk VariantFiltration -R {ref} -V $input -O $dir/$in_name"{params.gatk_VF_suf}.vcf" {params.gatk_VF_opt}
            input=$dir/$in_name"{params.gatk_VF_suf}.vcf"
            in_name=$in_name"{params.gatk_VF_suf}"
        fi
        gatk_SV_opt="{params.gatk_SV_opt}"
        if [[ ! -z $gatk_SV_opt ]]
        then
            echo -e "### filtering by gatk SelectVariants: {params.gatk_SV_opt}"
            gatk SelectVariants -R {ref} -V $input -O $dir/$in_name"{params.gatk_SV_suf}.vcf" {params.gatk_SV_opt}
            input=$dir/$in_name"{params.gatk_SV_suf}.vcf"
            in_name=$in_name"{params.gatk_SV_suf}"
        fi
        vcftools_opt="{params.vcftools_opt}"
        if [[ ! -z $vcftools_opt ]]
        then
            echo -e "### filtering by vcftools: {params.vcftools_opt}"
            cd $dir
            vcftools --vcf $in_name".vcf" --out $in_name"{params.vcftools_suf}" --recode {params.vcftools_opt}
            mv $in_name"{params.vcftools_suf}".recode.vcf $in_name"{params.vcftools_suf}".vcf
        fi
        """

rule concate_vcf:
    input: expand(os.path.join(output_dir, "{final_dir}/{chrom}{filtered}.vcf"), final_dir = final_dir, chrom = chrom, filtered = filtered)
    output:
        vcf_gz = os.path.join(output_dir, "{final_dir}/all_final{filtered}.vcf.gz"),
        vcf_tbi = os.path.join(output_dir, "{final_dir}/all_final{filtered}.vcf.gz.tbi"),
    shell:
        """
        out_vcf=$(echo {output.vcf_gz} | sed 's/\.gz//')
        bcftools concat -o $out_vcf {input}
        bgzip $out_vcf
        tabix -p vcf {output.vcf_gz}
        """


## stats
rule bcftools:
    input: os.path.join(output_dir, "{final_dir}/{chrom}{filtered}.vcf")
    output: os.path.join(output_dir, "stats/{final_dir}/{chrom}{filtered}.vcf")

rule test:
    shell:
        """
        echo "{x[1]}"
        """
# to make rule graph: snakemake gbp_variantcalling -np --rulegraph | dot -Tsvg > dag.svg
# to make job graph: snakemake gbp_variantcalling -np --dag | dot -Tsvg > dag_job.svg
