#!/bin/bash
#SBATCH --job-name=high_dim_job
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=100-102

####cd scratch/paper_iit_pt
module load r/4.4.0

Rscript --vanilla highd_example.R $SLURM_ARRAY_TASK_ID

