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
        vcf_gz = expand(os.path.join(output_dir, "{final_dir}/all_final{filtered}_singletons.vcf.gz"), final_dir = final_dir, filtered = filtered),
        vcf_tbi = expand(os.path.join(output_dir, "{final_dir}/all_final{filtered}_singletons.vcf.gz.tbi"), final_dir = final_dir, filtered = filtered),
        # expand(os.path.join(output_dir, "{final_dir}/{chrom}{filtered}.vcf"), final_dir = final_dir, chrom = chrom, filtered = filtered)

## pre-processing if neccessary

# rule fastp

## variant calling from fastq
rule index_ref:
    input: ref
    output: multiext(ref, ".bwt", ".pac", ".ann", ".amb", ".sa") # bwt. pac, ann, amb, sa
    conda: "conda.yaml"
    shell:
        """
        bwa index {ref}
        """


rule bwa_mem:
    input:
        fastq = os.path.join(fastq_dir, "{sample}"),
        ref_index = rules.index_ref.output
    output: temp(os.path.join(output_dir, "bwa_mem/{sample}.sam"))
    log: os.path.join(output_dir, "logs/snakemake/bwa_mem/{sample}.log")
    params: config["bwa_mem"]["params"]
    threads: config["bwa_mem"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        fastq=$(find {input.fastq} -regex ".*\(1\|2\).\(fq.gz\|fastq.gz\|fq\|fastq\)" | xargs -n1 | sort | xargs)
        echo "bwa mem -t {threads} {params} {ref} $fastq > {output}"
        bwa mem {ref} $fastq -t {threads} {params} -R $(./read_group.sh $fastq)> {output}
        """

rule picard_SortSam:
    input: rules.bwa_mem.output
    output:
        bam = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}.bam")),
        idx = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}.bai")),
    log: os.path.join(output_dir, "logs/snakemake/sorted_bam/picard_SortSam_{sample}.log")
    params: config["picard_SortSam"]["params"]
    conda: "conda.yaml"
    shell:
        """
        exec > >(tee {log}) 2>&1
        picard SortSam -I {input} -O {output.bam} -CREATE_INDEX true {params}
        """

rule samtools_view:
    input: rules.picard_SortSam.output.bam
    output:
        bam = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}_samtoolsView.bam")),
        idx = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}_samtoolsView.bam.bai")),
    log: os.path.join(output_dir, "logs/snakemake/sorted_bam/samtoolsView_{sample}.log")
    params: config["samtools_view"]["params"]
    threads: config["samtools_view"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        samtools view {params} -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        """

rule MarkDuplicates:
    input: rules.samtools_view.output.bam
    output:
        bam = os.path.join(output_dir, "sorted_bam/{sample}/{sample}_MarkDuplicate.bam"),
    log: os.path.join(output_dir, "logs/snakemake/sorted_bam/MarkDuplicate_{sample}.log")
    params: config["MarkDuplicates"]["params"]
    threads: config["MarkDuplicates"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        rm -rf {output}.parts
        ref_name=$(basename -s .fa {ref})
        ref_dict="$(dirname {ref})/$ref_name.dict"
        ref_fai="{ref}.fai"
        if [[ ! -f $ref_dict || ! -f $ref_fai ]]
        then
            picard CreateSequenceDictionary -REFERENCE {ref} -OUTPUT $ref_dict
            samtools faidx {ref}
        fi
        gatk MarkDuplicatesSpark -R {ref} -I {input} -O {output} --spark-runner LOCAL
        """

rule HaplotypeCaller:
    input: rules.MarkDuplicates.output
    output: temp(os.path.join(output_dir, "HaplotypeCaller/{sample}.g.vcf")),
    log: os.path.join(output_dir, "logs/snakemake/HaplotypeCaller/{sample}.log")
    params: config["HaplotypeCaller"]["params"]
    threads: config["HaplotypeCaller"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
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
    log: os.path.join(output_dir, "logs/snakemake/gendb/{chrom}.log")
    params: config["GenomicDBImport"]["params"]
    threads: config["GenomicDBImport"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        gatk GenomicsDBImport -V {GenomicDBImport_input} --genomicsdb-workspace-path {output} --intervals {wildcards.chrom} {params} --reader-threads {threads}
        """

rule GenotypeGVCFs :
    input: rules.GenomicDBImport.output
    output: os.path.join(output_dir, "vcf_by_chrom/{chrom}.vcf")
    log: os.path.join(output_dir, "logs/snakemake/vcf_by_chrom/{chrom}.log")
    params: config["CombineGVCFs"]["params"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        ref=$(realpath {ref})
        output=$(realpath {output})
        input=$(basename {input})
        cd $(dirname {input})
        gatk  GenotypeGVCFs  -R $ref -V gendb://$input -O $output {params} --intervals {wildcards.chrom}
        """



## variant filtration


rule filter_variants:
    input: rules.GenotypeGVCFs.output
    output: os.path.join(output_dir, "filtered_vcf_by_chrom/{chrom}{filtered}_singletons.vcf".format(chrom = "{chrom}", filtered = filtered))
    log: os.path.join(output_dir, "logs/snakemake/filtered_vcf_by_chrom/{chrom}.log")
    params:
        gatk_VF_opt = gatk_VF_opt,
        gatk_SV_opt = gatk_SV_opt,
        vcftools_opt = vcftools_opt,
        gatk_VF_suf = gatk_VF_suf,
        gatk_SV_suf = gatk_SV_suf,
        vcftools_suf = vcftools_suf,
        variant_count = "variant_count.csv",
    conda: "conda.yaml"
    shell:
        """
        exec > >(tee {log}) 2>&1
        init_dir=$PWD
        dir=$(dirname {output})
        stats_file=$dir/variant_count.csv
        if [[ ! -f $stats_file ]]
        then
            cp {params.variant_count} $dir/
            echo -e "Filter,None,{gatk_VF_opt},{gatk_SV_opt},{vcftools_opt},singletons_doubletons" >> $stats_file
            echo -e "Filename,None,{gatk_VF_suf},{gatk_SV_suf},{vcftools_suf},{vcftools_suf}_singletons" >> $stats_file
        fi
        raw_count=$(echo $(gatk CountVariants -V {input}) | sed 's/Tool returned://g')
        in_name=$(basename {input} | sed 's/.vcf//')
        input="{input}"
        gatk_VF_opt="{params.gatk_VF_opt}"
        if [[ ! -z $gatk_VF_opt ]]
        then
            echo -e "### filtering by gatk VariantFiltration: {params.gatk_VF_opt}"
            gatk VariantFiltration -R {ref} -V $input -O $dir/$in_name"{params.gatk_VF_suf}.vcf" {params.gatk_VF_opt}
            vf_count=$(echo $(gatk CountVariants -V $dir/$in_name"{params.gatk_VF_suf}.vcf") | sed 's/Tool returned://g')
            input=$dir/$in_name"{params.gatk_VF_suf}.vcf"
            in_name=$in_name"{params.gatk_VF_suf}"
        else
            vf_count=$raw_count
        fi
        gatk_SV_opt="{params.gatk_SV_opt}"
        if [[ ! -z $gatk_SV_opt ]]
        then
            echo -e "### filtering by gatk SelectVariants: {params.gatk_SV_opt}"
            gatk SelectVariants -R {ref} -V $input -O $dir/$in_name"{params.gatk_SV_suf}.vcf" {params.gatk_SV_opt}
            sv_count=$(echo $(gatk CountVariants -V $dir/$in_name"{params.gatk_SV_suf}.vcf") | sed 's/Tool returned://g')
            input=$dir/$in_name"{params.gatk_SV_suf}.vcf"
            in_name=$in_name"{params.gatk_SV_suf}"
        else
            sv_count=$vf_count
        fi
        vcftools_opt="{params.vcftools_opt}"
        if [[ ! -z $vcftools_opt ]]
        then
            echo -e "### filtering by vcftools: {params.vcftools_opt}"
            cd $dir
            vcftools --vcf $in_name".vcf" --out $in_name"{params.vcftools_suf}" --recode {params.vcftools_opt}
            mv $in_name"{params.vcftools_suf}".recode.vcf $in_name"{params.vcftools_suf}.vcf"
            vt_count=$(echo $(gatk CountVariants -V $in_name"{params.vcftools_suf}.vcf") | sed 's/Tool returned://g')
            cd $init_dir
        else
            vt_count=$sv_count
        fi
        cd $dir
        vcftools --singletons --vcf $in_name"{params.vcftools_suf}.vcf" --out "filter"
        if [[ $(cat filter.singletons | wc -l ) != 1 ]]
        then
            echo "### filtering by vcftools: singletons and doubletons"
            echo "singletons"
            vcftools --vcf $in_name"{params.vcftools_suf}.vcf" --out $in_name"{params.vcftools_suf}_singletons.vcf" --positions filter.singletons
            echo "filter singletons"
            sgl_count=$(echo $(gatk CountVariants -V $in_name"{params.vcftools_suf}_singletons.vcf") | sed 's/Tool returned://g')
        else
            cp $in_name"{params.vcftools_suf}.vcf" $in_name"{params.vcftools_suf}_singletons.vcf"
            sgl_count=$vt_count
        fi
        cd $init_dir
        echo -e "{wildcards.chrom},$raw_count,$vf_count,$sv_count,$vt_count,$sgl_count" >> $stats_file
        """

rule concate_vcf:
    input: expand(os.path.join(output_dir, "{final_dir}/{chrom}{filtered}_singletons.vcf"), final_dir = final_dir, chrom = chrom, filtered = filtered)
    output:
        vcf_gz = os.path.join(output_dir, "{final_dir}/all_final{filtered}.vcf.gz"),
        vcf_tbi = os.path.join(output_dir, "{final_dir}/all_final{filtered}.vcf.gz.tbi"),
    log: os.path.join(output_dir, "logs/snakemake/{final_dir}/concate_vcf_{filtered}.log")
    conda: "conda.yaml"
    shell:
        """
        exec > >(tee {log}) 2>&1
        out_vcf=$(echo {output.vcf_gz} | sed 's/\.gz//')
        bcftools concat -o $out_vcf {input}
        bgzip $out_vcf
        tabix -p vcf {output.vcf_gz}

        dir=$(dirname {output.vcf_gz})
        stats_file=$dir/variant_count.csv
        if [[ -f $stats_file ]]
        then
            append="Total"
            for i in 2 3 4 5 6
            do
                t=$(tail -n +4 $stats_file | cut -d',' -f$i | awk '{{Total=Total+$1}} END{{print Total}}')
                append=$append","$t
            done
            echo $append >> $stats_file
        fi
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
