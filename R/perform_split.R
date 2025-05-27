#' Perform a Split and Compute Objective (revised version)
#'
#' Attempts to perform a valid binary split on the given feature values, and evaluates the total objective value across the resulting child nodes.
#' If a proposed split produces one or more nodes with too few observations (less than \code{min.node.size}),
#' the function iteratively adjusts the split points by removing those causing small nodes, up to a maximum number of attempts.
#'
#' @param split.points A numeric vector of candidate split points (usually length 1).
#' @param xval A numeric vector of feature values to split on.
#' @param y A matrix of ICE values or any target values over which the objective function is computed.
#' @param min.node.size The minimum number of observations required in each resulting child node.
#' @param objective A function that returns a numeric value indicating the "loss" within a node, such as L2 loss or variance.
#'                  It must accept at least arguments \code{y}, \code{x}, and \code{requires.x}.
#' @param ... Additional arguments passed to the objective function.
#' @param max.attempts Maximum number of retries to find a valid split with all node sizes ≥ \code{min.node.size}.
#'
#' @return A single numeric value indicating the total objective value across valid child nodes, or \code{Inf} if no valid split could be found.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Mapping the samples to child nodes based on the split points using \code{findInterval()}.
#'   \item Checking whether each child node satisfies the minimum node size constraint.
#'   \item If not, it removes the offending split point(s) and retries, up to \code{max.attempts}.
#' }
#'
#' The splitting process avoids generating nodes with fewer than \code{min.node.size} observations. If no valid configuration is found,
#' the function returns \code{Inf} to indicate an invalid or poor split.
#'
#' @examples
#' \dontrun{
#'   perform_split(c(0.5), xval = runif(100), y = matrix(rnorm(1000), ncol = 10),
#'                 min.node.size = 10, objective = SS_L2)
#' }
#'
#' @seealso \code{\link{get_closest_point}}, \code{\link{find_best_binary_split}}
#'
perform_split <- function(split.points, xval, y, min.node.size, objective, ..., max.attempts = 5) {
  split.points = sort.int(split.points)
  split.points = get_closest_point(split.points, xval, min.node.size)

  attempt = 0
  while (attempt < max.attempts) {
    attempt <- attempt + 1
    # Put the samples into the interval according to the split points, get each sample belongs to the child node number,
    # because the original return value of `findInterval()` starts numbering from 0, so +1.
    node.number = findInterval(x = xval, split.points, rightmost.closed = TRUE) + 1
    node.size = tabulate(node.number) # Count the frequency of each unique value in the integer vector x, and return a vector of integers.

    if (min(node.size) >= min.node.size) {
      y.list = split(y, node.number)
      requires.x = formals(objective)[["requires.x"]]
      x.list = if (isTRUE(requires.x)) split(xval, node.number) else NULL

      res = vapply(seq_along(y.list), FUN = function(i) {
        objective(y = y.list[[i]], x = x.list[[i]], ...)
      }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)

      return(sum(res))
    }

    small.nodes = which(node.size < min.node.size) # Find out which intervals are “too small” and need to be handled.
    # If there are no more unreasonable nodes (small.nodes is empty), or the split points have been deleted (`split.points` is empty), then exit the while loop.
    if (length(small.nodes) == 0 || length(split.points) == 0) break
     # Map all unreasonable node numbers `small.nodes` to their corresponding left-hand split point indexes, and prevent 0-value indexes, de-duplicated and ready for deletion.
    drop.index = unique(pmax(small.nodes - 1, 1))
    # Ensure that the index of the split point to be deleted does not exceed the actual current number of splits
    drop.index = drop.index[drop.index <= length(split.points)]
    split.points = split.points[-drop.index]

    if (length(split.points) == 0) break
  }

  return(Inf)
}
