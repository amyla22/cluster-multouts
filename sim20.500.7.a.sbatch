#!/bin/bash

#SBATCH --partition=standard --time=120:00:00  --output=simout.%a.txt
#SBATCH --array=1-100
#SBATCH --mem=20GB

module load r/3.3.1/b1

# usage: Rscript <name of R script> <number of job from array> <input parameter file> <name of seed file>
Rscript SimRun.R $SLURM_ARRAY_TASK_ID SetParam.R Seedfile.txt
