#!/bin/bash
#SBATCH --job-name=hdim_A-IIT_p1k
#SBATCH --output=output_%a_%A.log
#SBATCH --error=error_%a_%A.log
#SBATCH --time=08:30:00
#SBATCH --mem=3GB
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=600-609

####cd scratch/paper_iit_pt
module load r/4.4.0

Rscript --vanilla codes/highd_example.R $SLURM_ARRAY_TASK_ID
