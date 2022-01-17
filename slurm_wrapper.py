#!/usr/bin/env python3

'''
wrap sbatch jobs based on Snakemake config file
'''
import os
import sys
from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
# config = sys.argv[1]

job_properties = read_job_properties(jobscript)


# config_properties = load_configfile(config)

rule = job_properties['rule']
job_name = '--job-name ' + rule

partition = '--partition long'

threads = job_properties['threads']
cpus_per_task = '--cpus-per-task ' + str(threads)

ntasks = '--ntasks 1'

log = str(job_properties['log'][0]) # os.path.join(output_dir, "logs/snakemake/bwa_mem/{sample}.log")
# logdir = os.path.dirname(os.path.dirname(os.path.dirname(log)))
# slurm_log = os.path.join(logdir, "slurm", os.path.basename(os.path.dirname(log)), os.path.basename(log))
# os.makedirs(os.path.dirname(slurm_log), exist_ok=True)

output = f'--output {log}_slurm_%j'
error = f'--error {log}_slurm_%j'

try:
    mem = job_properties['resources']['mem']
    mem = '--mem ' + str(mem)
except (IndexError, KeyError):
    mem = ''

cmdline = ' '.join(['sbatch --parsable ', job_name, partition, cpus_per_task, ntasks, mem, output, error, jobscript])


# jobscript.replace("\n", "\necho -e\"sbatch parameters:\n\"{}\"\"".format(sbatch), 1)
with open(jobscript, "r") as j:
    scripts = j.readlines()
scripts.insert(1, "echo -e \"# Submit command-line: \"{}\"\"\n".format(cmdline))
scripts.insert(2, "echo -e \"# Job running on node: $SLURM_JOB_NODELIST\"\n")
scripts.insert(3, "echo -e \"\n\"")
with open(jobscript, "w") as j:
    j.writelines(scripts)

os.system(cmdline)
# print(cmdline)
# sbatch --job-name {cluster.job-name} --partition {cluster.partition} --account {cluster.account} --cpus-per-task {cluster.cpus-per-task} --output {cluster.output} --error {cluster.error} snakejob.sh
