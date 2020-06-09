library('PAIRADISE')

args <- commandArgs(trailingOnly=TRUE)

input_file_name <- args[1]
number_of_threads_str <- args[2]
output_file_name <- args[3]

data_frame <- read.table(file=input_file_name, sep="\t", header=TRUE)

number_of_threads <- as.integer(number_of_threads_str)

# pairadise has an error for an input data set with only 1 event.
# This is because an intermediate result in the calculation ends up as
# a 1-dimensional vector if there is only 1 input event, but a
# 2-dimensional table is expected.
#
# Avoid the 1 event error by adding an extra event to the input before
# running pairadise and then removing that extra event before
# writing to the output file.

num_input_rows <- nrow(data_frame)
# For an empty input, just add empty PValue and FDR columns.
if (num_input_rows == 0) {
  data_frame["PValue"] <- vector()
  data_frame["FDR"] <- vector()
} else {
  if (num_input_rows == 1) {
    filler_df <- data.frame("FillerID", "1", "0", "0", "1", "20", "10")
    names(filler_df) <- c("ID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2",
                          "SJC_SAMPLE_2", "IncFormLen", "SkipFormLen")
    data_frame <- rbind(data_frame, filler_df)
  }
  results <- PAIRADISE::pairadise(data_frame,
                                  numCluster=number_of_threads)

  # pairadise may filter out some events. Fill in NA for those rows.
  num_rows <- nrow(data_frame)
  num_results <- length(results$raw.pvalues)
  fdr_values <- p.adjust(results$raw.pvalues, method="BH")
  p_value_column <- vector(mode="numeric", length=num_rows)
  fdr_column <- vector(mode="numeric", length=num_rows)
  results_i <- 1
  for (df_i in 1:num_rows) {
    if ((results_i <= num_results)
        && (data_frame$ID[df_i] == results$exonID[results_i])) {
      p_value_column[df_i] <- results$raw.pvalues[results_i]
      fdr_column[df_i] <- fdr_values[results_i]
      results_i <- results_i + 1
    } else {
      p_value_column[df_i] <- NA
      fdr_column[df_i] <- NA
    }
  }

  data_frame["PValue"] <- p_value_column
  data_frame["FDR"] <- fdr_column

  if (num_input_rows == 1) {
    data_frame <- data_frame[1, ]
  }
}

write.table(data_frame,
            file=output_file_name,
            sep="\t",
            quote=FALSE,
            row.names=FALSE)
