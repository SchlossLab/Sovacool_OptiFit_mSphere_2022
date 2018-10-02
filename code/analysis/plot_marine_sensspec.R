library(tidyverse)
make_plot <- function() {
  sensspec <- read.table(file = "data/marine/marine.sensspec.final", header = TRUE)
  sensspec <- sensspec[, 4:ncol(sensspec)]
  
  sensspec_long <- gather(sensspec, key = "Score", value = "Value", c(1:11, 13, 14))
  
  #Plot of mcc values against the % of data used as reference
  mcc <- ggplot(data = sensspec, aes(x = refp, y = mcc)) +
    geom_jitter(width = 0.1, height = 0) +
    labs(title = "OptiFIt de novo clustering using partial human dataset") +
    ylab("MCC") +
    xlab("Amount of data used as reference") +
    #scale_y_continuous(breaks = seq(0, 1, by = .1))+
    #coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 100, by = 10))
  
  ggsave("results/figures/marine/marine_mcc.png", plot = mcc, width = 6.5, height = 4.5, unit="in")
  
  #Plot mcc values against various other scores from the confusion matrix
  mcc_v_cm <- ggplot(data=sensspec_long, aes(x = Value, y = mcc)) +
    geom_point() +
    ylab("MCC") +
    xlab("Amount of data used as reference") +
    facet_wrap(~ Score, scales="free")
  
  ggsave("results/figures/marine/marine_mcc_v_cm.png", plot=mcc_v_cm, width = 10, height = 7, unit="in")
}