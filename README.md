## My SNP calling pipeline from fastq files to a single analysis-ready VCF

### following GATK best practice for Germline short variant discovery
[<img src="https://drive.google.com/uc?id=1HKtzOeobgOVjCXEUE0-5378ocBz6Age7">](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

### my pipeline
`snakemake gbp_variantcalling --use-conda`

<img src="./dag_job.svg">
