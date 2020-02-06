# Adapted from: https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017/blob/master/code/soil.batch
library(tidyverse)

make_files_file <- function(metadata_filename, output_filename){
  samples <- read_tsv(metadata_filename) %>%
    select(Run_s) %>%
    mutate(read_1 = paste0(Run_s, '_1.fastq.gz'),
           read_2 = paste0(Run_s, '_1.fastq.gz'))
  write_tsv(samples, output_filename, col_names = FALSE)
}

make_files_file("data/soil/soil.metadata", "data/soil/soil.files")
