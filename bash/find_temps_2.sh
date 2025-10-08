#!/bin/bash
#SBATCH --job-name=find_t_PT-IIT_G
#SBATCH --output=output_%a_%A.log
#SBATCH --error=error_%a_%A.log
#SBATCH --time=48:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=168,170

module load r/4.5.0

Rscript --vanilla codes/find_temp_manual.R $SLURM_ARRAY_TASK_ID
