#Uses a connections file and an accnos file to see how many highly connected sequences are in the reference
#Usage: check_connections.R infile
args = commandArgs(trailingOnly = TRUE)
infile = args[1]

check_connections <- function(infile) {
  #Get all sequences and take the top 10% most connected ones
  connections <- read.table(file = infile, header = TRUE)
  connections <- connections[1:floor(nrow(connections)/10), ]
  connections$v1 <- as.character(connections$v1)
  
  #Sequences that were chosen for the sample
  accnos <- read.table(file = "data/marine/sample.accnos", header = FALSE)
  accnos <- as.character(accnos$V1)
  
  #We want to know how many connected sequences are in the reference, which means
  #we want to know how many are NOT in the sample
  ref_seqs <- sum(!(connections$v1 %in% accnos))
  
  #Cat here prevents R from printing [1] in front of the number
  cat(ref_seqs)
}

check_connections(infile)