require(ggplot2)
require(readr)
require(dplyr)

#Makes MCC plots from the provided sensspec file
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
#Make the plots
make_plot <- function(filename) {
  #Read in accumulated sensspec data created by optifit_multi.sh
  sensspec <- readr::read_tsv(file = filename) %>%
    dplyr::select(c(4:ncol(.))) #First 3 columns don't contain any useful information
  
  combo_mcc <- ggplot2::ggplot(data=sensspec, ggplot2::aes(x = refp, y = mcc, color = type)) +
    ggplot2::geom_jitter(width = 1, height = 0) +
    ggplot2::labs(title = "OptiFit clustering using partial marine dataset") +
    ggplot2::ylab("MCC") +
    ggplot2::xlab("Amount of data used as reference") +
    #ggplot2::scale_y_continuous(breaks = seq(0, 1, by = .1))+
    #ggplot2::coord_cartesian(ylim = c(.75, 1)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 10))
  
  ggplot2::ggsave(paste("results/figures/marine/", basename(filename), ".mcc.png", sep = ""),
                  plot = combo_mcc, width = 6.5, height = 4.5, unit = "in")
  
  combo_mcc_full <- ggplot2::ggplot(data=sensspec, ggplot2::aes(x = refp, y = mcc, color = type)) +
    ggplot2::geom_jitter(width = 1, height = 0) +
    ggplot2::labs(title = "OptiFit clustering using partial marine dataset") +
    ggplot2::ylab("MCC") +
    ggplot2::xlab("Amount of data used as reference") +
    #ggplot2::scale_y_continuous(breaks = seq(0, 1, by = .1))+
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 10))
  
  ggplot2::ggsave(paste("results/figures/marine/", basename(filename), ".mcc.full.png", sep = ""),
                  plot = combo_mcc_full, width = 6.5, height = 4.5, unit = "in")
}

make_plot(filename)