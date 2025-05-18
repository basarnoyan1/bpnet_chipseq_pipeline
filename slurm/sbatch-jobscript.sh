#!/bin/bash
#SBATCH --job-name=smk
#SBATCH --output=logs/slurm/job-%j.out
#SBATCH --error=logs/slurm/job-%j.err
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=mid
#SBATCH --gpus-per-node=1
#SBATCH --mem=256G
#SBATCH --qos=users
#SBATCH --account=users
#SBATCH --time=4:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bnoyan21@ku.edu.tr

# Load modules
module load anaconda/2024.02
source activate basar_test_py39
module load apptainer/1.3.4
module load cudnn/8.9.5/cuda-11.x
module load cuda/11.8.0

# Navigate to pipeline directory
export PATH="/kuacc/users/bnoyan21/hpc_run/sratoolkit.3.1.1-ubuntu64/bin:$PATH"
cd /kuacc/users/bnoyan21/hpc_run/Tasks/Task130525

# Run the Snakemake job
{exec_job}
