#!/bin/bash
#SBATCH --job-name=high_dim_job
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=88-92

####cd scratch/paper_iit_pt
module load r/4.4.0

Rscript --vanilla highd_example.R $SLURM_ARRAY_TASK_ID

