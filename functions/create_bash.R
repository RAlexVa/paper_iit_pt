create_sh_files <- function(numbers,
                            prefix = "high",
                            hours=48,
                            cpus=20,
                            num_sims=100,
                            output_dir = ".",
                            name_file="default") {
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
#SBATCH --mem=5GB
#SBATCH --cpus-per-task=CPUS
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=1-NUMSIMS

module load r/4.5.0

FIXED_INPUT=PLACEHOLDER

Rscript --vanilla codes/highd_seeded.R $FIXED_INPUT $SLURM_ARRAY_TASK_ID
'
    # Replace the placeholder with the actual number
    content <- gsub("PLACEHOLDER", numbers[i], content)
    content <- gsub("HOURS", hours, content)
    content <- gsub("CPUS", cpus, content)
    content <- gsub("NUMSIMS", num_sims, content)
    content <- gsub("NAMEFILE", name_file, content)
    
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
numbers <- c(846:853)

create_sh_files(numbers,
                prefix = "high",
                hours=100,
                cpus=10,
                num_sims=20,
                output_dir = "bash",
                name_file="1k_6modes")

numbers <- c(854:877)
create_sh_files(numbers,
                prefix = "higher",
                hours=5,
                cpus=20,
                num_sims=50,
                output_dir = "bash",
                name_file="3-5-7k_6modes")
numbers <- c(878:883)
create_sh_files(numbers,
                prefix = "repli",
                hours=20,
                cpus=20,
                num_sims=50,
                output_dir = "bash",
                name_file="3-5k_6modes")