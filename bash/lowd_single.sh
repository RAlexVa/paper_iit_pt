#!/bin/bash
#SBATCH --job-name=lowd_single
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=12
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=FAIL,BEGIN,END

PARAM="209"


module load r/4.4.0

Rscript --vanilla codes/lowd_example.R $PARAM

