#' Map Split Points to Closest Valid Data Values
#'
#' For each candidate split point, finds the closest actual data value from a reduced set of
#' well-spaced points in \code{xval} that ensure balanced node sizes.
#'
#' @param split.points A numeric vector of proposed split thresholds.
#' @param xval A numeric vector of feature values (should be sorted).
#' @param min.node.size Minimum number of observations required on both sides of a split.
#'
#' @return A sorted numeric vector of split points mapped to nearest valid values.
#'
#' @details
#' This step avoids duplicate or unstable splits by ensuring that:
#' \itemize{
#'   \item Returned values exist in the actual data.
#'   \item Each value respects node size constraints.
#'   \item Values are reused at most once.
#' }
#'
#' @seealso \code{\link{perform_split}}
#'
get_closest_point = function(split.points, xval, min.node.size = 10) {
  xval = sort.int(xval)
  # Ensure `min.node.size` between points (is not guaranteed if many duplicated values exist)
  chunk.ind = seq.int(min.node.size + 1, length(xval) - min.node.size, by = min.node.size)
  xadj = unique(xval[chunk.ind]) # unique(quantile(xval, prob = chunk.ind/length(xval), type = 1))
  # xval = xval[-c(1, length(xval))]
  split.adj = numeric(length(split.points))
  # Iterate over all split points, replacing them with the closest `xadj` value.
  for (i in seq_along(split.adj)) {
    d = xadj - split.points[i]
    ind.closest = which.min(abs(d))
    split.adj[i] = xadj[ind.closest]
    xadj = xadj[-ind.closest] # remove already chosen value
  }

  return(sort.int(split.adj))
}
