#!/bin/bash
#SBATCH --job-name=A-IIT_2k_20t
#SBATCH --output=output_%a_%A.log
#SBATCH --error=error_%a_%A.log
#SBATCH --time=36:30:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=670-679

####cd scratch/paper_iit_pt
module load r/4.4.0

Rscript --vanilla codes/highd_example.R $SLURM_ARRAY_TASK_ID
