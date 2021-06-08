#' No Duplicate Dataset
#'
#' A dataset containing \code{730} simulated records from \code{3} files with
#' no duplicate records in each file.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{records}{A \code{data.frame} with the records, containing \code{7}
#'   fields, from all three files, in the format used for input to
#'   \code{\link{create_comparison_data}}.}
#'   \item{file_sizes}{The size of each file.}
#'   \item{IDs}{The true partition of the records, represented as an
#'   \code{integer}  vector of arbitrary labels of length
#'   \code{sum(file_sizes)}.}
#' }
#' @source Extracted from the datasets used in the simulation study of the
#' paper. The datasets were generated using code from Peter Christen's group
#' \url{https://dmm.anu.edu.au/geco/index.php}.
#' @examples
#' data(no_dup_data)
#'
#' # There are 500 entities represented in the records
#' length(unique(no_dup_data$IDs))
"no_dup_data"
