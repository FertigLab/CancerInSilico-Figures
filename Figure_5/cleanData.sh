#!/bin/bash -l

#SBATCH
#SBATCH --job-name=fig5_clean
#SBATCH --time=0:30:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=tsherma4@jhu.edu

if [ "$1" != "" ]; then
    DIR=$1
else
    echo "no input directory"
    exit 1
fi

Rscript cleanData.R $DIR
