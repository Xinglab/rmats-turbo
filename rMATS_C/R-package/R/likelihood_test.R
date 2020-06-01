## Examples of C code

#'
#' @export
#' @useDynLib rMATS


lt = function(inputf, cutoff, nthread=1) {
    row_num <- length(readLines(file(inputf, open='r')))-1
    gc()
    res <- .C('r_likelihood_test', as.character(inputf), as.double(cutoff), as.integer(nthread), double(row_num), package = 'rMATS')
    res[[4L]]
}
