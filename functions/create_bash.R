create_sh_files <- function(numbers,
                            prefix = "high",
                            hours=48,
                            cpus=20,
                            sim_ini=1,
                            sim_fin=100,
                            memory_gb=16,
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
#SBATCH --mem=MEMORYGB
#SBATCH --cpus-per-task=CPUS
#SBATCH --mail-user=alexander.valencia@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --array=SIMINI-SIMFIN

module load r/4.5.0

FIXED_INPUT=PLACEHOLDER

Rscript --vanilla codes/highd_seeded.R $FIXED_INPUT $SLURM_ARRAY_TASK_ID
'
    # Replace the placeholder with the actual number
    content <- gsub("PLACEHOLDER", numbers[i], content)
    content <- gsub("HOURS", hours, content)
    content <- gsub("CPUS", cpus, content)
    content <- gsub("SIMINI", sim_ini, content)
    content <- gsub("SIMFIN", sim_fin, content)
    content <- gsub("NAMEFILE", name_file, content)
    content <- gsub("MEMORY", memory_gb, content)
    
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
numbers <- c(1058:1061)
create_sh_files(numbers,
                prefix = "spa1",
                hours=100,
                cpus=20,
                sim_ini=1,
                sim_fin=20,
                memory_gb=25,
                output_dir = "bash",
                name_file="spaced_500")
numbers <- c(1062:1065)
create_sh_files(numbers,
                prefix = "spa2",
                hours=100,
                cpus=20,
                sim_ini=21,
                sim_fin=40,
                memory_gb=25,
                output_dir = "bash",
                name_file="spaced_500")




numbers <- c(1052:1057)
create_sh_files(numbers,
                prefix = "high",
                hours=20,
                cpus=20,
                sim_ini=1,
                sim_fin=25,
                memory_gb=2,
                output_dir = "bash",
                name_file="1k_13temp_diftemp")

numbers <- c(1012:1023)
create_sh_files(numbers,
                prefix = "less_rf",
                hours=8,
                cpus=20,
                sim_ini=1,
                sim_fin=25,
                memory_gb=5,
                output_dir = "bash",
                name_file="6mod_13t_less_rf")

numbers <- c(978:1007)
create_sh_files(numbers,
                prefix = "more_rf",
                hours=100,
                cpus=20,
                sim_ini=1,
                sim_fin=25,
                memory_gb=16,
                output_dir = "bash",
                name_file="6mod_13t_less_rf")