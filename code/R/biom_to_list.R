# original source: https://github.com/SchlossLab/Schloss_Cluster_PeerJ_2015/blob/master/code/biom_to_list.R
# Modified Oct. 2020 by KLS

library(dplyr)
library(jsonlite)

biom_to_list <- function(biom_file_name, list_file_name, label = 0.03) {
  biom_file <- scan(biom_file_name, what = "", sep = "\n", quiet = TRUE)
  biom_json <- fromJSON(biom_file)

  n_otus <- biom_json$shape[1]
  n_seqs <- biom_json$shape[2]
  seq_names <- biom_json$columns[["id"]]
  otu_names <- biom_json$rows[["id"]]

  biom_json$data <- biom_json$data + 1
  otu_assignments <- aggregate(seq_names[biom_json$data[, 2]],
    by = list(biom_json$data[, 1]),
    function(x) {
      paste(x, collapse = ",")
    }
  )$x

  # storing this in memory as a giant string is a bad idea,
  # but I only have to run it once per vsearch clustering job,
  # so I don't really care.
  list_data <- paste(paste(c("label", "numOTUs", otu_names),
    collapse = "\t"
  ),
  paste(c(label, n_otus, otu_assignments),
    collapse = "\t"
  ),
  sep = "\n"
  )

  write(list_data, list_file_name)
}

biom_to_list(snakemake@input[["biom"]], snakemake@output[["list"]])
