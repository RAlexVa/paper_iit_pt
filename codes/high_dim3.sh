#!/bin/bash
#SBATCH --job-name=high_dim_job
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=28:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

####cd scratch/paper_iit_pt
module load r/4.4.0

PARAM="99"
Rscript --vanilla codes/highd_example.R $PARAM

