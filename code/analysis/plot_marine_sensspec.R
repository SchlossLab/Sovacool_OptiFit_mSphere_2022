#Makes MCC plots from the provided sensspec file
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]

#Make the plots
make_plot <- function(filename) {
  #Read in accumulated sensspec data created by optifit_multi.sh
  sensspec <- read.table(file = filename, header = TRUE)
  sensspec <- sensspec[, 4:ncol(sensspec)]
  sensspec <- tidyr::gather(sensspec, key = "mcc_type", value = "mcc", c(12, 15, 16))
  
  combo_mcc <- ggplot2::ggplot(data=sensspec, ggplot2::aes(x = refp, y = mcc, color = mcc_type)) +
    ggplot2::geom_jitter(width = 1, height = 0) +
    ggplot2::labs(title = "OptiFit clustering using partial marine dataset") +
    ggplot2::ylab("MCC") +
    ggplot2::xlab("Amount of data used as reference") +
    #ggplot2::scale_y_continuous(breaks = seq(0, 1, by = .1))+
    ggplot2::coord_cartesian(ylim = c(.75, 1)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 10))
  
  ggplot2::ggsave(paste("results/figures/marine/", basename(filename), ".mcc.png", sep = ""),
                  plot = combo_mcc, width = 6.5, height = 4.5, unit = "in")
  
  sensspec_long <- tidyr::gather(sensspec, key = "Score", value = "Value", c(1:11, 13, 14))
  
  #Plot mcc values against various other scores from the confusion matrix
  mcc_v_cm <- ggplot2::ggplot(data=sensspec_long, ggplot2::aes(x = Value, y = mcc)) +
    ggplot2::geom_point() +
    ggplot2::ylab("MCC") +
    ggplot2::xlab("Amount of data used as reference") +
    ggplot2::facet_wrap(~ Score, scales="free")
  
  ggplot2::ggsave(paste("results/figures/marine/", basename(filename), ".mcc_v_cm.png", sep = ""),
                  plot=mcc_v_cm, width = 10, height = 7, unit="in")
}

make_plot(filename)