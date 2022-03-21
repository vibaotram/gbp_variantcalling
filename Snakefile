import os
from Bio import SeqIO
import re
import glob
import math

configfile: "config.yaml"

## input
fastq_dir = config["fastq_dir"]
samples,fq, = glob_wildcards(os.path.join(fastq_dir, "{sample}/{fq, .*.(fq|fq.gz|fastq|fastq.gz)}"))
uniq_samples = set(samples)
sample=[]
for i in uniq_samples:
    sample.append(os.path.basename(i))
fastq_files, = glob_wildcards(os.path.join(fastq_dir, "{fastq_files, .*/.*.(fq|fq.gz|fastq|fastq.gz)}"))

ref = config["ref"]
chrom = []
# fasta_sequences = SeqIO.parse(open(ref),'fasta')
# for fasta in fasta_sequences:
#     seqname = fasta.id
#     if not re.match("Contig", seqname):
#         chrom.append(seqname)

records = list(SeqIO.parse(open(ref),'fasta'))
idlist = list()
for r in records:
    idlist.append(r.id)
chrom = []
for s in idlist:
    if not re.match("Contig", s):
        chrom.append(s)
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
# rule fastqc:
#     input: os.path.join(output_dir, "fastqc/fastqc.html")

rule fastqc:
    input: expand(os.path.join(fastq_dir, "{fastq_files}"), fastq_files = fastq_files)
    output: os.path.join(output_dir, "fastqc/fastqc_final.html")
    log: os.path.join(output_dir, "logs/fastqc/fastqc.log")
    params:
        opt = config["fastqc"]["params"]
    threads: config["fastqc"]["threads"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        outdir=$(dirname {output})
        mkdir -p $outdir
        fastqc -o $outdir -t {threads} -f fastq {input}
        multiqc -o $outdir -n fastqc $outdir/*_fastqc.zip
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
    log: os.path.join(output_dir, "logs/bwa_mem/index_ref.log")
    shadow: "full"
    conda: "conda.yaml"
    shell:
        """
        exec > >(tee {log}) 2>&1
        bwa index {ref}
        """

def input_bwa_mem(wildcards):
    # d,fq, = glob_wildcards(os.path.join(fastq_dir, "{d}/{wildcards.sample}/{fq, .*.(fq|fq.gz|fastq|fastq.gz)}"))
    fq = glob.glob("{fastq_dir}/**/{sample}/*.f*".format(fastq_dir = fastq_dir, sample = wildcards.sample), recursive=True)
    return fq


rule bwa_mem:
    input:
        fastq = input_bwa_mem,
        # fastq = os.path.join(fastq_dir, "{sample}"),
        ref_index = rules.index_ref.output
    output: temp(os.path.join(output_dir, "bwa_mem/{sample}.sam"))
    shadow: "full"
    log: os.path.join(output_dir, "logs/bwa_mem/{sample}.log")
    params: config["bwa_mem"]["params"]
    threads: config["bwa_mem"]["threads"]
    resources:
        mem = config["bwa_mem"]["mem"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        # fastq=$(find {input.fastq} -regex ".*\(1\|2\).\(fq.gz\|fastq.gz\|fq\|fastq\)" | xargs -n1 | sort | xargs)
        echo "bwa mem -t {threads} {params} {ref} {input.fastq} > {output}"
        bwa mem {ref} {input.fastq} -t {threads} {params} -R $(./read_group.sh {input.fastq}) > {output}
        """

rule picard_SortSam:
    input: rules.bwa_mem.output
    output:
        bam = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}_SortSam.bam")),
        bai = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}_SortSam.bai")),
    shadow: "full"
    log: os.path.join(output_dir, "logs/sorted_bam/picard_SortSam_{sample}.log")
    params: config["picard_SortSam"]["params"]
    resources:
        mem = config["picard_SortSam"]["mem"]
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
        bai = temp(os.path.join(output_dir, "sorted_bam/{sample}/{sample}_samtoolsView.bam.bai")),
    shadow: "full"
    log: os.path.join(output_dir, "logs/sorted_bam/samtoolsView_{sample}.log")
    params: config["samtools_view"]["params"]
    threads: config["samtools_view"]["threads"]
    resources:
        mem = config["samtools_view"]["mem"]
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
        bai = os.path.join(output_dir, "sorted_bam/{sample}/{sample}_MarkDuplicate.bam.bai"),
        sbi = os.path.join(output_dir, "sorted_bam/{sample}/{sample}_MarkDuplicate.bam.sbi"),
    shadow: "full"
    log: os.path.join(output_dir, "logs/sorted_bam/MarkDuplicate_{sample}.log")
    params: config["MarkDuplicates"]["params"]
    threads: config["MarkDuplicates"]["threads"]
    resources:
        mem = config["MarkDuplicates"]["mem"]
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
        gatk MarkDuplicatesSpark -R {ref} -I {input} -O {output.bam} --spark-runner LOCAL
        """

rule HaplotypeCaller:
    input: rules.MarkDuplicates.output.bam
    output:
        gvcf = os.path.join(output_dir, "HaplotypeCaller/{sample}.g.vcf"),
        idx = os.path.join(output_dir, "HaplotypeCaller/{sample}.g.vcf.idx"),
    shadow: "full"
    log: os.path.join(output_dir, "logs/HaplotypeCaller/{sample}.log")
    params: config["HaplotypeCaller"]["params"]
    threads: config["HaplotypeCaller"]["threads"]
    resources:
        mem = config["HaplotypeCaller"]["mem"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        if [[ ! -f {ref}.fai ]]
        then
            samtools faidx {ref} --fai-idx {ref}.fai
        fi
        gatk HaplotypeCaller -R {ref} -I {input} -O {output.gvcf} --native-pair-hmm-threads {threads} -ERC GVCF {params}
        """

GenomicDBImport_input = ' -V '.join(expand(rules.HaplotypeCaller.output, sample = uniq_samples))



def GenomicDBImport_intervals(wildcards):
    chrid = wildcards.chrom
    chridx = idlist.index(chrid)
    size = 1e7
    chrlen = len(records[chridx])
    chr_itv = math.floor(chrlen/size)
    if chrlen % size < 5e5:
        chr_itv = chr_itv
    else:
        chr_itv += 1
    GenomicDBImport_intervals = ''
    for i in range(0, chr_itv):
        f = int((i * size) + 1)
        if i < chr_itv -1:
            l = int((i + 1) * size)
        else:
            l = chrlen
        s = '-L {chrid}:{f}-{l} '.format(chrid=chrid, f=f, l=l)
        GenomicDBImport_intervals += s
    return GenomicDBImport_intervals



rule GenomicsDBImport:
    input: expand(rules.HaplotypeCaller.output.gvcf, sample = set(sample))
    output: temp(directory(os.path.join(output_dir, "gendb/{chrom}")))
    shadow: "full"
    log: os.path.join(output_dir, "logs/gendb/{chrom}.log")
    params:
        intervals = GenomicDBImport_intervals,
        ext_params = config["GenomicsDBImport"]["params"]
    # threads: config["GenomicsDBImport"]["threads"]
    resources:
        mem = config["GenomicsDBImport"]["mem"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        gatk GenomicsDBImport -V {GenomicDBImport_input} --genomicsdb-workspace-path {output} --overwrite-existing-genomicsdb-workspace -imr OVERLAPPING_ONLY {params.intervals}{params.ext_params}
        """

# rule CombineGVCFs:
#     input: expand(rules.HaplotypeCaller.output.gvcf, sample = sample)
#     output:
#         gvcf = os.path.join(output_dir, "vcf_by_chrom/combined.g.vcf"),
#         idx = os.path.join(output_dir, "vcf_by_chrom/combined.g.vcf.idx"),
#     shadow: "full"
#     log: os.path.join(output_dir, "logs/vcf_by_chrom/CombineGVCFs.log")
#     params: config["CombineGVCFs"]["params"]
#     resources:
#         mem = config["CombineGVCFs"]["mem"]
#     conda: "conda.yaml"
#     # singularity: singularity_img
#     shell:
#         """
#         exec > >(tee {log}) 2>&1
#         input=""
#         for i in {input}; do input=$input' -V '$i; done
#         gatk CombineGVCFs -R {ref} $input -O {output.gvcf} {params}
#         """

# rule GenotypeGVCFs :
#     input: rules.CombineGVCFs.output.gvcf
#     output:
#         gvcf = os.path.join(output_dir, "vcf_by_chrom/{chrom}.vcf"),
#         idx = os.path.join(output_dir, "vcf_by_chrom/{chrom}.vcf.idx"),
#     shadow: "full"
#     log: os.path.join(output_dir, "logs/vcf_by_chrom/GenotypeGVCFs_{chrom}.log")
#     params: config["GenotypeGVCFs"]["params"]
#     resources:
#         mem = config["GenotypeGVCFs"]["mem"]
#     conda: "conda.yaml"
#     # singularity: singularity_img
#     shell:
#         """
#         exec > >(tee {log}) 2>&1
#         gatk GenotypeGVCFs -R {ref} -V {input} -O {output.gvcf} {params} --intervals {wildcards.chrom}
#         """

rule GenotypeGVCFs :
    input: rules.GenomicsDBImport.output
    output:
        gvcf = os.path.join(output_dir, "vcf_by_chrom/{chrom}.vcf"),
        idx = os.path.join(output_dir, "vcf_by_chrom/{chrom}.vcf.idx"),
    shadow: "full"
    log: os.path.join(output_dir, "logs/vcf_by_chrom/GenotypeGVCFs_{chrom}.log")
    params: config["GenotypeGVCFs"]["params"]
    resources:
        mem = config["GenotypeGVCFs"]["mem"]
    conda: "conda.yaml"
    # singularity: singularity_img
    shell:
        """
        exec > >(tee {log}) 2>&1
        cd $(dirname {input})
        db=$(basename {input})
        gatk GenotypeGVCFs -R {ref} -V gendb://$db -O {output.gvcf} {params} --intervals {wildcards.chrom}
        """

## variant filtration


rule filter_variants:
    input: rules.GenotypeGVCFs.output.gvcf
    output: os.path.join(output_dir, "filtered_vcf_by_chrom/{chrom}{filtered}_singletons.vcf".format(chrom = "{chrom}", filtered = filtered))
    shadow: "full"
    log: os.path.join(output_dir, "logs/filtered_vcf_by_chrom/{chrom}.log")
    params:
        gatk_VF_opt = gatk_VF_opt,
        gatk_SV_opt = gatk_SV_opt,
        vcftools_opt = vcftools_opt,
        gatk_VF_suf = gatk_VF_suf,
        gatk_SV_suf = gatk_SV_suf,
        vcftools_suf = vcftools_suf,
        variant_count = "variant_count.csv",
    resources:
        mem = config["filtration_mem"]
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
            cd $dir
            vcftools --vcf $in_name"{params.gatk_VF_suf}.vcf" --out $in_name"{params.gatk_VF_suf}" --remove-filtered-all --recode --recode-INFO-all
            mv $in_name"{params.gatk_VF_suf}.recode.vcf" $in_name"{params.gatk_VF_suf}.vcf"
            cd $init_dir
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
            vcftools --vcf $in_name".vcf" --out $in_name"{params.vcftools_suf}" --recode --recode-INFO-all {params.vcftools_opt}
            mv $in_name"{params.vcftools_suf}.recode.vcf" $in_name"{params.vcftools_suf}.vcf"
            vt_count=$(echo $(gatk CountVariants -V $in_name"{params.vcftools_suf}.vcf") | sed 's/Tool returned://g')
            cd $init_dir
        else
            vt_count=$sv_count
        fi
        cd $dir
        vcftools --singletons --vcf $in_name"{params.vcftools_suf}.vcf" --out {wildcards.chrom}
        if [[ $(cat "{wildcards.chrom}.singletons" | wc -l ) != 1 ]]
        then
            echo "### filtering by vcftools: singletons and doubletons"
            vcftools --vcf $in_name"{params.vcftools_suf}.vcf" --out $in_name"{params.vcftools_suf}_singletons" --positions "{wildcards.chrom}.singletons" --recode --recode-INFO-all
            mv $in_name"{params.vcftools_suf}_singletons.recode.vcf" $in_name"{params.vcftools_suf}_singletons.vcf"
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
    shadow: "full"
    log: os.path.join(output_dir, "logs/{final_dir}/concate_vcf{filtered}.log")
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

rule test:
    shell:
        """
        echo "{samples}"
        """
# to make rule graph: snakemake gbp_variantcalling -np --rulegraph | dot -Tsvg > dag.svg
# to make job graph: snakemake gbp_variantcalling -np --dag | dot -Tsvg > dag_job.svg
