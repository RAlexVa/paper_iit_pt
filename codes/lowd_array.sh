#!/bin/bash
#SBATCH --job-name=lowd_final_res
#SBATCH --output=output_%a_%A.log
#SBATCH --error=error_%a_%A.log
#SBATCH --time=10:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=FAIL ###BEGIN,END,

#SBATCH --array=201-216,219-222

####cd scratch/paper_iit_pt
### we can use percentage j to get the job ID
module load r/4.4.0

Rscript --vanilla codes/lowd_example.R $SLURM_ARRAY_TASK_ID

