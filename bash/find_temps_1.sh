#!/bin/bash
#SBATCH --job-name=find_t_PT-IIT_1k
#SBATCH --output=output_%a_%A.log
#SBATCH --error=error_%a_%A.log
#SBATCH --time=20:00:00
#SBATCH --mem=5GB
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=89-91,83-85

module load r/4.4.0

Rscript --vanilla codes/find_temp_manual.R $SLURM_ARRAY_TASK_ID
