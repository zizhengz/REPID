#' Find the Best Binary Split for a Single Feature
#'
#' Evaluates a set of candidate split points for a single feature using a specified objective function,
#' and returns the split point with the lowest loss.
#'
#' @param xval A numeric vector of feature values.
#' @param y A matrix or data.frame of responses (e.g. ICE curves).
#' @param n.splits Number of splits to consider (must be \code{1} for binary trees).
#' @param min.node.size Minimum number of observations per resulting node.
#' @param objective A loss function to be minimized, e.g., \code{SS_L2}. It must support arguments \code{y}, \code{x}, and optionally \code{requires.x}.
#' @param ... Extra arguments forwarded to \code{objective} or downstream utilities.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{split.points}}{The best-performing split point (scalar).}
#'   \item{\code{objective.value}}{The minimized loss achieved at that split.}
#' }
#'
#' @details
#' Internally:
#' \itemize{
#'   \item Candidates are generated via \code{\link{generate_split_candidates}}.
#'   \item Each candidate is evaluated using \code{\link{perform_split}}.
#' }
#'
#' @seealso \code{\link{perform_split}}, \code{\link{generate_split_candidates}}
#'
find_best_binary_split <- function(xval, y, n.splits = 1, min.node.size = 10, objective, ...) {
  # `n.splits` must be equal to 1, otherwise an error is reported -> binary tree
  assert_choice(n.splits, choices = 1)

  # Use different split candidates to perform split
  q = generate_split_candidates(xval, n.quantiles = 100, min.node.size = min.node.size)
  splits = vapply(q, FUN = function(i) {
    perform_split(i, xval = xval, y = y, min.node.size = min.node.size,
                  objective = objective, ...)
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  # Select the split point yielding the minimal objective
  best = which.min(splits)

  return(list(split.points = q[best], objective.value = splits[best]))
}
