require(ggplot2)
require(readr)
require(dplyr)

#Usage: call script from snakemake workflow
#Filename must be the full path from working directory to the sensspec file
#dataset is the name of the dataset being used (human, mice, marine soil)
#Makes MCC plots from the provided sensspec fit

make_plot <- function(dataset, input_filename, output_filenames) {
  title <- paste("OptiFit clustering using", dataset, "dataset", sep = " ")
  
  #Read in accumulated sensspec data
  sensspec <- readr::read_tsv(file = input_filename) %>%
    dplyr::select(c(4:ncol(.))) #First 3 columns don't contain any useful information
  
  combo_mcc <- ggplot2::ggplot(data=sensspec, ggplot2::aes(x = reference_fraction, y = mcc, color = type)) +
    ggplot2::geom_jitter(width = 1, height = 0) +
    ggplot2::labs(title = title) +
    ggplot2::ylab("MCC") +
    ggplot2::xlab("Amount of data used as reference (%)") +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 10))
  
  ggplot2::ggsave(output_filenames[["combo_mcc"]],
                  plot = combo_mcc, width = 6.5, height = 4.5, unit = "in")
  
  combo_mcc_full <- ggplot2::ggplot(data=sensspec, ggplot2::aes(x = reference_fraction, y = mcc, color = type)) +
    ggplot2::geom_jitter(width = 1, height = 0) +
    ggplot2::labs(title = title) +
    ggplot2::ylab("MCC") +
    ggplot2::xlab("Amount of data used as reference (%)") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 10))
  
  ggplot2::ggsave(output_filenames[["mcc_full"]],
                  plot = combo_mcc_full, width = 6.5, height = 4.5, unit = "in")
  
  iters <- ggplot2::ggplot(data=dplyr::filter(sensspec, reference_fraction == .5 & (type == "SAMP" | type == "SAMP_O_NOREF")),
                           ggplot2::aes(x = reference_fraction, y = mcc, color = type)) +
    ggplot2::geom_jitter(width = .1, height = 0, size = 3) +
    ggplot2::labs(title = title) +
    ggplot2::ylab("MCC") +
    ggplot2::xlab("Cut") +
    ggplot2::scale_x_continuous(breaks = seq(1, 10, by = 1))
  
  ggplot2::ggsave(output_filenames[["iters"]],
                  plot = iters, width = 6.5, height = 4.5, unit = "in")
}

dataset <- snakemake@params[["dataset"]]
input_filename <- snakemake@input[[1]]
output_filenames <- snakemake@output
make_plot(dataset, input_filename, output_filenames)