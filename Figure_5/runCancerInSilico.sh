#!/bin/bash -l

#SBATCH
#SBATCH --job-name=runCIS
#SBATCH --time=3:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=tsherma4@jhu.edu

if [ "$1" != "" ]; then
    FILE=$1
else
    echo missing file name
    exit 1
fi

if [ "$SLURM_ARRAY_TASK_ID" != "" ]; then
    ARRAY_NUM=$SLURM_ARRAY_TASK_ID
else
    ARRAY_NUM=1
fi

time Rscript $FILE $ARRAY_NUM
