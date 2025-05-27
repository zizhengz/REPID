#' Optimize Feature Split for ICE Tree Node
#'
#' Evaluates all candidate split points across features and returns the one with the best objective value.
#' Used to determine the best way to split a node when building interpretable trees for ICE curves.
#'
#' @param Y A matrix or \code{data.table} of response values (typically centered ICE curves).
#' @param X A \code{data.frame} or \code{data.table} of features used for splitting.
#' @param n.splits Number of split points to consider (must be 1 for binary splits).
#' @param min.node.size Minimum number of observations required in each resulting node after split.
#' @param optimizer A function to find the best binary split for a feature, e.g. \code{find_best_binary_split()}.
#' @param objective A function that computes the loss to minimize (e.g., L2 or area loss on ICE).
#' @param ... Additional arguments passed to the optimizer and objective function.
#'
#' @return A \code{data.table} with one row per feature, containing:
#' \describe{
#'   \item{\code{feature}}{The name of the splitting feature.}
#'   \item{\code{objective.value}}{The minimum loss achieved using this feature.}
#'   \item{\code{runtime}}{Elapsed time in seconds to evaluate splits on the feature.}
#'   \item{\code{split.points}}{A list-column of selected split points (typically of length 1).}
#'   \item{\code{best.split}}{A logical flag indicating whether this feature achieved the best (lowest) objective.}
#' }
#'
#' @details
#' This function loops through each column in \code{X}, applies the given \code{optimizer} function
#' to find the best split for that feature, and collects results in a ranked \code{data.table}.
#' Only binary splits (i.e., \code{n.splits = 1}) are currently supported.
#'
#' @seealso \code{\link{find_best_binary_split}}, \code{\link{compute_tree}}
#'
split_parent_node = function(Y, X, n.splits = 1, min.node.size = 10, optimizer,
                             objective, ...) {

  assert_data_frame(X)
  assert_integerish(n.splits)
  assert_integerish(min.node.size)
  assert_function(objective, args = c("y", "x", "requires.x"))
  assert_function(optimizer, args = c("xval", "y"))

  # Find best split points per feature
  opt.feature = lapply(X, function(feat) {
    t1 = proc.time()
    res = optimizer(x = feat, y = Y, n.splits = n.splits, min.node.size = min.node.size,
                    objective = objective, ...)
    t2 = proc.time()
    res$runtime = (t2 - t1)[[3]]
    return(res)
  })

  result = data.table::rbindlist(lapply(opt.feature, as.data.frame), idcol = "feature")
  result = result[, .(split.points = list(split.points)), by = c("feature", "objective.value", "runtime"), with = TRUE]
  result$best.split = result$objective.value == min(result$objective.value)
  return(result)
}
