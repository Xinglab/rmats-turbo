library('Darts')

args <- base::commandArgs(trailingOnly=TRUE)

input_file_name <- args[1]
output_file_name <- args[2]
number_of_threads_str <- args[3]
cutoff_str <- args[4]
has_replicates_str <- args[5]

number_of_threads <- base::as.integer(number_of_threads_str)
cutoff <- base::as.numeric(cutoff_str)
has_replicates <- has_replicates_str == 'true'

num_input_lines <- base::length(base::readLines(input_file_name, n=2))
## If the input is empty or only has a header then create an empty output file
if (num_input_lines < 2) {
    out_connection <- base::file(output_file_name, 'w')
    base::close(out_connection)
} else {
  if (has_replicates) {
    Darts::Darts_replicate(input_file_name, output_file_name, C=cutoff,
                           thread=number_of_threads, estim_groupVar_prior=FALSE)
  } else {
    Darts::Darts(input_file_name, output_file_name, C=cutoff,
                 thread=number_of_threads)
  }
}
