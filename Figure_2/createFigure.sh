#!/bin/bash -l

#SBATCH
#SBATCH --job-name=fig2
#SBATCH --time=0:30:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=18
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=tsherma4@jhu.edu

Rscript createData.R 18
Rscript createFigure.R
