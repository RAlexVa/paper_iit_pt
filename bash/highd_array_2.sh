#!/bin/bash
#SBATCH --job-name=PT-IIT_1k_7t
#SBATCH --output=output_%a_%A.log
#SBATCH --error=error_%a_%A.log
#SBATCH --time=16:30:00
#SBATCH --mem=3GB
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=660-669

####cd scratch/paper_iit_pt
module load r/4.4.0

Rscript --vanilla codes/highd_example.R $SLURM_ARRAY_TASK_ID
