#' Small Duplicate Dataset
#'
#' A dataset containing \code{53} simulated records from \code{3} files with
#' no duplicate records in each file, subset from \code{\link{dup_data}}.
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
#' @references Serge Aleshin-Guendel & Mauricio Sadinle (2022). Multifile Partitioning for Record Linkage and Duplicate Detection. \emph{Journal of the
#' American Statistical Association}. [\doi{https://doi.org/10.1080/01621459.2021.2013242}][\href{https://arxiv.org/abs/2110.03839}{arXiv}]
#'
#' @examples
#' data(dup_data_small)
#'
#' # There are 53 entities represented in the records
#' length(unique(dup_data_small$IDs))
"dup_data_small"
