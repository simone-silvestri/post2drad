#!/bin/bash -l
#SBATCH --time=120:00:00
##SBATCH --time=01:00:00
#SBATCH -N 1 
#SBATCH --ntasks-per-node=1
##SBATCH --partition=short
#SBATCH --output="out_s"
#SBATCH --mem=54000

MV2_ENABLE_AFFINITY=0

scp simone@int2-bb.cartesius.surfsara.nl:/archive/simone/backup-scratch/develop_cooled/DATA/* data/
