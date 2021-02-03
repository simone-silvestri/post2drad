#!/bin/bash -l
#
#SBATCH --nodes=96
#SBATCH --ntasks=1152
#SBATCH --ntasks-per-node=12
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=gpu
#SBATCH --time=24:00:00
#SBATCH --output="out_s"
#SBATCH --account=s1002
#SBATCH --mem=62GB

# tasks: $SLURM_NTASKS 
# tasks-per-node: $SLURM_NTASKS_PER_NODE 
# cpus-per-task: $SLURM_CPUS_PER_TASK
srun --cpu_bind=rank ./postOutput

