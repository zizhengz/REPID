#' Generate Valid Candidate Split Points
#'
#' Produces a set of valid candidate split points from a numeric feature vector, ensuring that each candidate allows
#' at least \code{min.node.size} observations on both sides after splitting.
#'
#' @param xval A numeric vector of feature values to generate splits from.
#' @param n.quantiles Optional integer. If provided, candidate splits are based on quantile thresholds of the central region of \code{xval}.
#' @param min.node.size Minimum number of observations required per child node.
#'
#' @return A numeric vector of adjusted split point candidates that do not coincide with raw \code{xval} values,
#'         and ensure sufficient node size if used for splitting.
#'
#' @details
#' This function supports two modes:
#' \itemize{
#'   \item If \code{n.quantiles} is specified, splits are drawn as quantiles of the inner region of \code{xval}.
#'   \item Otherwise, evenly spaced valid midpoints are returned directly.
#' }
#' In both cases, split points are post-processed using \code{\link{adjust_split_point}} to avoid alignment with actual data values,
#' which improves split stability and reproducibility.
#'
#' @seealso \code{\link{adjust_split_point}}, \code{\link{find_best_binary_split}}
#'
generate_split_candidates <- function(xval, n.quantiles = NULL, min.node.size = 10) {
  # Check that `min.node.size` is an integer or an “approximate integer” and that its value does not exceed the lower bound of `(length(xval)-1)/2`.
  # This prevents insufficient data on either side of the split.
  assert_integerish(min.node.size, lower = 1, upper = floor((length(xval) - 1)/2))

  xval = sort.int(xval)
  # Ensure that only points in the middle are selected that satisfy the split condition, i.e. there are `min.node.size` observations on both left and right.
  chunk.ind = seq.int(min.node.size + 1, length(xval) - min.node.size, by = min.node.size)
  # Extract these candidate cutoffs from the sorted `xval`
  xadj = xval[chunk.ind]

  if (!is.null(n.quantiles)) { # Whether to generate split points based on quantiles
    # To speedup we use only quantile values as possible split points
    # qprobs = seq(1/n.quantiles, (n.quantiles - 1)/n.quantiles, by = 1/n.quantiles)
    qprobs = seq(0, 1, by = 1/n.quantiles)
    q = unique(quantile(xadj, qprobs, type = 1)) # `type = 1` -> empirical quantiles; `unique()` ensures that split points are not duplicated.
  } else {
    q = unique(xadj)
  }
   # The final result is a set of processed candidate split points that do not overlap with the values in the `xval`.
  q = adjust_split_point(q, xval)
  return(q)
}
