#!/bin/bash
#SBATCH --job-name=smk-{rule}
#SBATCH --output=logs/{rule}.{wildcards}.%j.out
#SBATCH --error=logs/{rule}.{wildcards}.%j.err
#SBATCH --cpus-per-task={threads}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH --mem=256G
#SBATCH --qos=users
#SBATCH --account=users
#SBATCH --time=10:00:00
#SBATCH --output=test-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnoyan21@ku.edu.tr


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