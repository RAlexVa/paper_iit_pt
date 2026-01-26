create_sh_files <- function(numbers,
                            prefix = "high",
                            hours=48,
                            cpus=20,
                            sim_ini=1,
                            sim_fin=100,
                            memory_gb=16,
                            output_dir = ".",
                            name_file="default",
                            name_code="highd_seeded.R") {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (i in seq_along(numbers)) {
    filename <- file.path(output_dir, paste0(prefix, "_", i, ".sh"))
    
    content <- '#!/bin/bash
#SBATCH --job-name=NAMEFILE
#SBATCH --output=output_%a_%A.log
#SBATCH --error=error_%a_%A.log
#SBATCH --time=HOURS:00:00
#SBATCH --mem=MEMORYGB
#SBATCH --cpus-per-task=CPUS
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=SIMINI-SIMFIN

module load r/4.5.0

FIXED_INPUT=PLACEHOLDER

Rscript --vanilla codes/CODE_NAME $FIXED_INPUT $SLURM_ARRAY_TASK_ID
'
    # Replace the placeholder with the actual number
    content <- gsub("PLACEHOLDER", numbers[i], content)
    content <- gsub("HOURS", hours, content)
    content <- gsub("CPUS", cpus, content)
    content <- gsub("SIMINI", sim_ini, content)
    content <- gsub("SIMFIN", sim_fin, content)
    content <- gsub("NAMEFILE", name_file, content)
    content <- gsub("MEMORY", memory_gb, content)
    content <- gsub("CODE_NAME", name_code, content)
    # # Write the file with UNIX line endings
    # writeLines(content, filename, sep = "\n")
    # # Make the file executable (Unix/Linux systems)
    # Sys.chmod(filename, mode = "755")
    
    conx = file(filename, open="wb")
    write(content, conx)
    close(conx)
    
    
    cat("Created:", filename, "\n")
  }
}

# Usage with your example vector
numbers <- c(700:703)
create_sh_files(numbers,
                prefix = "lowd",
                hours=3,
                cpus=1,
                sim_ini=1,
                sim_fin=100,
                memory_gb=6,
                output_dir = "bash",
                name_file="lowdim",
                name_code="lowd_seeded.R")

numbers <- c(1100:1103)
create_sh_files(numbers,
                prefix = "highd1",
                hours=10,
                cpus=20,
                sim_ini=1,
                sim_fin=100,
                memory_gb=6,
                output_dir = "bash",
                name_file="highdim_1k",
                name_code="highd_seeded.R")

numbers <- c(1104:1115)
create_sh_files(numbers,
                prefix = "highd",
                hours=8,
                cpus=20,
                sim_ini=1,
                sim_fin=100,
                memory_gb=6,
                output_dir = "bash",
                name_file="highdim",
                name_code="highd_seeded.R")


numbers <- c(704:705)
create_sh_files(numbers,
                prefix = "lowd2",
                hours=3,
                cpus=1,
                sim_ini=1,
                sim_fin=100,
                memory_gb=6,
                output_dir = "bash",
                name_file="lowdim",
                name_code="lowd_seeded.R")
