#Look at how "connected" each sequence is within the dataset
#A connection is defined as a similarity with another sequence of >97% (also <3% difference)

#Read in entire distance table
#This file contains pairwise distances between sequences
human_dists <- read.table("Data\\Human\\Processed\\human.0_4.01.sm.dist", header = FALSE,
                          col.names = c("Seq1", "Seq2", "Dist"))

#Remove distances >3%
human_dists <- human_dists[human_dists$Dist < 0.03, ]

#All pairwise distances only occur once, and a given seq could occur in either the first
#or second column
#Stack both columns together to get a count of all unique occurences of a sequence
human_connect_counts <- c(as.character(human_dists$Seq1), as.character(human_dists$Seq2))

#Count occurences of each sequence
human_connect_counts <- aggregate(data.frame(count = human_connect_counts),
                                  list(value = human_connect_counts),
                                  length)
human_connect_counts <- human_connect_counts[order(human_connect_counts$count, decreasing = TRUE), ]

#Plot showing raw counts of sequences with a given number of connections
human_connect_plot <- ggplot(data = human_connect_counts, aes(x = count)) +
  geom_histogram(bins = 500) +
  ylab("Number of sequences") +
  xlab("Number of pairwise connections")

human_connect_plot

#Cumulative distribution showing approximately how many sequences are "highly connected"
human_ecdf_plot <- ggplot(data = human_connect_counts, aes(x = count)) +
  stat_ecdf() +
  scale_y_continuous(breaks = seq(0, 1, by = .1)) +
  ylab("Number of sequences") +
  xlab("Number of pairwise connections")

human_ecdf_plot