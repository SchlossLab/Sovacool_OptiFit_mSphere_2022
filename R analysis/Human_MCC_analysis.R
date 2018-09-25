#First exploration of output from the OptiFit algorithm in Mothur 1.41.0
#Checking how variation in the percent of data used to create the reference
#affects the final MCC

library(dplyr) # For %>% operator
library(ggplot2)

#Read a Mothur sensspec file and returns it as a dataframe
#Also adds a column listing the percent of original data used as reference
read_sensspec <- function(filename, percent) {
  temp_df <- read.table(file = filename, sep = "\t", header = TRUE)
  temp_df$Percent <- percent
  temp_df
}

#Create filenames with full filepath
human_filenames <- dir("Data\\Human\\OptiFit Output") %>%
  lapply(function(x) paste("Data\\Human\\OptiFit Output\\", x, sep = ""))

#Percent of data used as a reference by optifit denovo clustering algorithm
#Our data was run from 5% to 95% at increments of 5%
#This data is not included in Mothur's output so must be provided manually
percents <- seq(95, 5, by = -5) 

#Read each file into a dataframe with read_sensspec
human_results <- lapply(1:19, function(x) read_sensspec(filename = human_filenames[[x]],
                                                        percent = percents[[x]])) 

human_results <- do.call("rbind", human_results) #Collapse list into a single df

#Plot MCC value of each run against the % of data used as the reference when de novo clustering
human_mccs <- ggplot(data = human_results, aes(x = Percent, y = mcc)) +
  geom_jitter(width = 0.1, height = 0) +
  labs(title = "OptiFIt de novo clustering using partial human dataset") +
  ylab("MCC") +
  xlab("Amount of data used as reference")

human_mccs

#Convert data to long form with all confusion matrix measures having their own row
human_cm_data <- reshape2::melt(human_results, id.vars = c("iter", "label", "cutoff", "Percent", "mcc"),
                      variable.name = "measure", value.name = "value")

human_cm_plot <- ggplot(data = human_cm_data, aes(x = value, y = mcc, color = Percent)) +
  geom_point() +
  facet_wrap(~ measure, scales = "free")

human_cm_plot