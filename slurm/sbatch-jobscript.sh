#!/bin/bash
#SBATCH --job-name=smk-{rule}
#SBATCH --output=logs/{rule}.{wildcards}.%j.out
#SBATCH --error=logs/{rule}.{wildcards}.%j.err
#SBATCH --time={resources.time}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={resources.mem_mb}
{resources.gpus and f"#SBATCH --gpus={resources.gpus}"}

#Basar's custom cluster settings
module load anaconda/2024.02
source activate basar_test_py39
module load apptainer/1.3.4
module load cudnn/8.9.5/cuda-11.x
module load cuda/11.8.0
cd /kuacc/users/bnoyan21/hpc_run/
export PATH=$PATH:$PWD/sratoolkit.3.1.1-ubuntu64/bin
cd fastq_files

{exec_job}