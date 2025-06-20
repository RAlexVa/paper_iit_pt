#!/bin/bash
#SBATCH --job-name=hdim_PT-IIT_p2k
#SBATCH --output=output_%a_%A.log
#SBATCH --error=error_%a_%A.log
#SBATCH --time=25:30:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=630-649

####cd scratch/paper_iit_pt
module load r/4.4.0

Rscript --vanilla codes/highd_example.R $SLURM_ARRAY_TASK_ID