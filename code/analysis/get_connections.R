library(dplyr)
#Reads a distance file and finds the number of pairwise distances that
#every sequence in the file has

#Returns a sorted table with the most connected sequences at the top
get_connections <- function() {
  distances <- read.table(file = "data/marine/marine.200.dist", header = FALSE)
  
  #Each sequence can appear in either column of an entry in the dist table,
  #so to get total number of connections per seq stack both columns together
  seqs <- c(as.character(distances$V1), as.character(distances$V2)) %>%
    factor() %>%
    data.frame()
  names(seqs) <- c("v1")
  seqs <- count(seqs, v1) %>%
    arrange(desc(n))
  seqs$v1 <- as.character(seqs$v1)
  
  write.table(seqs, file = "data/marine/marine.200.connections", row.names = FALSE, quote = FALSE)
}