# Load required packages
library(R6)
library(data.table)
library(assertthat)
library(ggplot2)
library(tidyr)
library(scales)
library(libcoin)
library(partykit)


#------------------------------------------------------------------------------------------------------------
# Functions for computing








#' Compute L2 Loss for ICE Curves
#'
#' Calculates the sum of squared deviations between individual ICE curves and their mean (PDP).
#'
#' @param y A matrix of ICE curve values.
#' @param x Not used in L2. Included for compatibility.
#' @param requires.x Logical, whether the objective requires x
#' @return Numeric total L2 loss.
#'
#' @export
SS_L2 <- function(y, x, requires.x = FALSE, ...) {
  ypred = colMeans(as.matrix(y))
  sum(t((t(y) - ypred)^2))
}






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
#' @export
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
#' @export
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








#' Slightly Adjust Split Points to Avoid Data Overlap
#'
#' Nudges candidate split points slightly so they fall between observed values of \code{xval},
#' avoiding exact alignment with actual data points. This helps prevent ambiguous node assignments during splitting.
#'
#' @param split.points A numeric vector of proposed split points, possibly overlapping with values in \code{xval}.
#' @param xval A numeric vector of observed feature values (must be sorted or sortable).
#'
#' @return A numeric vector of adjusted split points that lie strictly between values in \code{xval}, avoiding overlap.
#'
#'  @details
#' This function is typically used to:
#' \itemize{
#'   \item Prevent exact equality between split points and feature values, which can cause undefined or inconsistent behavior during splitting.
#'   \item Ensure the resulting intervals are strictly non-overlapping and well-defined.
#' }
#'
#' The adjustment works as follows:
#' \enumerate{
#'   \item If a split point matches an exact value in \code{xval}, nudge it rightward by half the gap to the next unique value.
#'   \item If not matched, shift it slightly by half the minimal difference in \code{xval}.
#'   \item Correct boundary cases to ensure all points fall between valid inner values (not near the edges).
#' }
#'
#' @seealso \code{\link{generate_split_candidates}}
#'
#' @export
adjust_split_point <- function(split.points, xval) {
  # `split.points`: original split points, may overlap with values in xval
  q = split.points
  # Take out unique values and sort them
  x.unique = sort.int(unique(xval))
  # Find the locations in `x.unique` of the split points that coincide with the unique values of `xval`
  ind = which(x.unique %in% q)
  # Remove the index of the last position (so that `ind+1` in the subsequent step cannot be out of bounds)
  ind = ind[ind < length(x.unique)]
  # If some values in the split points `q` are not in `x.unique` (not found by `ind`), use `min(diff(...)) /2` as the minimum interval for fine-tuning;
  # Otherwise, for each element in `q`, find the distance between it and the next unique value, taking half.
  eps = if (length(ind) != length(q)) min(diff(x.unique))/2 else (x.unique[ind + 1] - x.unique[ind])/2
  # Use `eps` to make all the split points nudge a little to the right and fall between the two observations.
  q = q + eps

  # Fix boundary problem: if one element in `q` is too far to the left or too far to the right after fine-tuning, adjust to the midpoint of the two neighboring unique values.
  q[q < x.unique[2]] = mean(x.unique[1:2])
  q[q > x.unique[length(x.unique) - 1]] = mean(x.unique[(length(x.unique) - 1):length(x.unique)])
  return(q)
}









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
#' @export
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










# Performs a single split and measures the objective
# -------
# perform_split <- function(split.points, xval, y, min.node.size, objective, ...) {
#   split.points = sort.int(split.points)
#   split.points = get_closest_point(split.points, xval, min.node.size)
#
#   # assign intervalnr. according to split points
#   node.number = findInterval(x = xval, split.points, rightmost.closed = TRUE) + 1 #把样本根据切分点划入区间,得到每个样本所属的子节点编号
#   # compute size of each child node
#   node.size = tabulate(node.number) # tabulate: 统计整数向量x中每个唯一值的频数,返回一个整数向量
#   # if minimum node size is violated, return Inf
#   # TODO: instead of returning Inf try to avoid that this happens by fixing split points
#   if (min(node.size) < min.node.size) return(Inf)
#   y.list = split(y, node.number) # compute objective in each interval and sum it up(按照切分结果,把所有y值分配到对应的子节点中,每个子节点接下来就可以独立计算误差了)
#   requires.x = formals(objective)[["requires.x"]] # 检查用户传入的objective函数是否需要使用x(特征值)作为输入变量之一
#   # x.list only needed if this is used in the objective
#   x.list = if (isTRUE(requires.x)) split(xval, node.number) else NULL
#   res = vapply(seq_along(y.list), FUN = function(i) {
#     objective(y = y.list[[i]], x = x.list[[i]], ...)
#   }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
#   sum(res)
# }
# ----


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
#' @export
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
#' @export
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





#' R6 Class: Node
#'
#' The \code{Node} class defines a recursive data structure used to build an interpretable decision tree over ICE curves.
#' Each node represents a subset of the data and stores information about its depth, objective value,
#' split decision, and its children nodes.
#'
#' @section Fields:
#' \describe{
#'   \item{\code{id}}{Integer identifier for the node. Not globally unique.}
#'   \item{\code{depth}}{Integer. Depth of the node in the tree (root = 1).}
#'   \item{\code{subset.idx}}{Integer vector of row indices from the dataset belonging to this node.}
#'   \item{\code{objective.value}}{Objective function value for the node.}
#'   \item{\code{objective.value.parent}}{Parent node's objective value (used for relative improvement).}
#'   \item{\code{split.feature}}{Feature name used for splitting (if split occurred).}
#'   \item{\code{split.value}}{Threshold used to split this feature.}
#'   \item{\code{children}}{A named list with \code{left.child} and \code{right.child}, each a \code{Node} or \code{NULL}.}
#'   \item{\code{stop.criterion.met}}{Logical. Whether the node met a stopping criterion and cannot be split.}
#'   \item{\code{improvement.met}}{Logical. Whether the relative improvement was below threshold.}
#'   \item{\code{intImp}}{Numeric. Improvement measure used for early stopping.}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{initialize(...)}}{Creates a new \code{Node} object with specified index subset and optional parent linkage.}
#'   \item{\code{computeSplit(X, Y, objective, impr.par, optimizer, min.split)}}{Attempts to split this node using the best split point across features.}
#'   \item{\code{computeChildren(X, Y, objective)}}{Constructs child \code{Node} objects for left and right partitions if split occurred.}
#' }
#'
#' @seealso \code{\link{compute_tree}}, \code{\link{split_parent_node}}
#'
#' @export
Node <- R6Class("Node", list(
  id = NULL,
  depth = NULL,# on which depth is the node
  subset.idx = NULL,# ids of the instances of data that are in this node
  objective.value = NULL, # objective value in a node
  objective.value.parent = NULL,
  id.parent = NULL,
  child.type = NULL,# left or right type

  # Split information (if splitting has already taken place)
  split.feature = NULL,
  split.value = NULL,

  # Append the children of this node
  children = list(),

  stop.criterion.met = FALSE,
  improvement.met = NULL,
  intImp = NULL,

  initialize = function(id, depth = NULL, subset.idx, id.parent = NULL, child.type = NULL,
                        objective.value.parent = NULL, objective.value = NULL, improvement.met, intImp) {

    assert_numeric(id, len = 1)
    assert_numeric(depth, len = 1, null.ok = TRUE)

    assert_numeric(subset.idx, min.len = 1)
    assert_numeric(id.parent, len = 1, null.ok = TRUE)
    assert_character(child.type, null.ok = TRUE)

    self$id = id
    self$depth = depth
    self$subset.idx = subset.idx
    self$id.parent = id.parent
    self$child.type = child.type
    self$intImp = intImp
    self$objective.value.parent = objective.value.parent
    self$objective.value = objective.value
    self$stop.criterion.met = FALSE
    self$improvement.met = improvement.met
  },

  computeSplit = function(X, Y, objective, impr.par, optimizer, min.split = 10) {
    if (length(self$subset.idx) < min.split | self$improvement.met == TRUE) {
      self$stop.criterion.met = TRUE
    } else {
      self$objective.value.parent = objective(y = Y, x = X)
      self$objective.value = objective(y = Y[self$subset.idx, ], x = X[self$subset.idx, ])
      tryCatch({
        split = split_parent_node(Y = Y[self$subset.idx, ], X = X[self$subset.idx, ],
                                  objective = objective, optimizer = optimizer, min.node.size = min.split)

        if (is.null(self$intImp)) self$intImp = 0
        intImp = (self$objective.value - split$objective.value[split$best.split][1]) / self$objective.value.parent

        if ((self$intImp == 0 && intImp < impr.par) || (self$intImp != 0 && intImp < self$intImp * impr.par)) {
          self$improvement.met = TRUE
        } else {
          self$split.feature = split$feature[split$best.split][1]
          self$split.value = unlist(split$split.points[split$best.split])[1]
          self$intImp = intImp
          self$objective.value.parent = objective(y = Y[self$subset.idx, ], x = X[self$subset.idx, ])
          self$objective.value = split$objective.value[split$best.split][1]
        }
      },
      error = function(e) {
        self$stop.criterion.met = TRUE
      })
    }
  },

  computeChildren = function(X, Y, objective) {
    if (self$stop.criterion.met | self$improvement.met) {
      # no further split is performed
      self$children = list("left.child" = NULL, "right.child" = NULL)
    } else {
      if(is.null(self$split.feature))
        stop("Please compute the split first via computeSplit().")

      idx.left = which(X[self$subset.idx, self$split.feature, with = FALSE] <= self$split.value)
      idx.right = which(X[self$subset.idx, self$split.feature, with = FALSE] > self$split.value)

      idx.left = self$subset.idx[idx.left]
      idx.right = self$subset.idx[idx.right]

      if (length(idx.left) == 0) idx.left = 0
      if (length(idx.right) == 0) idx.right = 0

      #obj.left = objective(y = Y[idx.left, ], x = X[idx.left, ])
      #obj.right = objective(y = Y[idx.right, ], x = X[idx.right, ])
      #obj.parent = objective(y = Y[self$subset.idx, ], x = X[self$subset.idx, ])

      left.child = Node$new(id = 1, depth = self$depth + 1, subset.idx = idx.left,
                            id.parent = self$id, child.type = "<=", improvement.met = self$improvement.met,
                            intImp = self$intImp)
      right.child = Node$new(id = 2, depth = self$depth + 1, subset.idx = idx.right,
                             id.parent = self$id, child.type = ">", improvement.met = self$improvement.met,
                             intImp = self$intImp)
      self$children = list("left.child" = left.child, "right.child" = right.child)
    }
  }
))










#' Build the ICE-Based Decision Tree
#'
#' Constructs a binary tree over ICE curves using recursive node splits.
#'
#' @param effect An object of class \code{FeatureEffect} from the \pkg{iml} package. Must contain ICE data.
#' @param testdata A \code{data.frame} or \code{data.table} of predictor variables used to generate ICE curves.
#' @param objective Character. One of \code{"SS_L1"}, \code{"SS_L2"}, \code{"SS_area"}, or \code{"SS_sd"}.
#'                 Specifies the objective function to use for splitting.
#' @param n.split Integer. Maximum depth (number of levels) of the resulting tree.
#' @param impr.par Numeric. Minimum relative improvement required to continue splitting a node.
#' @param min.split Integer. Minimum number of samples required in a node to be eligible for splitting.
#' @return A list-of-lists representing the ICE regression tree. Each element corresponds to a level of the tree,
#'         containing \code{Node} objects or \code{NULL} (for terminated branches).
#'
#' @details
#' The function initializes with a root node containing all samples and recursively applies
#' \code{Node$computeSplit()} and \code{Node$computeChildren()} to expand the tree.
#'
#' Objective options:
#' \describe{
#'   \item{\code{"SS_L1"}}{Minimum squared distance from ICE curves to the best-fitting curve.}
#'   \item{\code{"SS_L2"}}{Total squared deviation from the PDP (mean of ICE).}
#'   \item{\code{"SS_area"}}{Variance of area-under-curve per sample.}
#'   \item{\code{"SS_sd"}}{Total variance of model predictions (used for standard deviation trees).}
#' }
#'
#' @seealso \code{\link{Node}}, \code{\link{compute_data_for_ice_splitting}}, \code{\link{find_best_binary_split}}
#'
#' @export
compute_tree = function(effect, testdata, objective = "SS_L2", n.split, impr.par = 0.05, min.split = 10) {

  if (objective == "SS_L1") {
    split.objective = function(y, x, requires.x = FALSE, ...) {
      require(Rfast)
      ypred = colMeans(as.matrix(y)) # column means → PDP
      min(t((t(y) - ypred)^2)) # Represented by an ICE curve that most closely approximates average behavior
    }
    input.data = compute_data_for_ice_splitting(effect, testdata = testdata)
  }
  else if (objective == "SS_L2") {
    split.objective = function(y, x, requires.x = FALSE, ...) {
      ypred = colMeans(as.matrix(y)) # PDP
      sum(t((t(y) - ypred)^2)) # Total squared deviation from the PDP
    }
    input.data = compute_data_for_ice_splitting(effect, testdata = testdata)
  }
  else if (objective == "SS_area") {
    split.objective = function(y, x, requires.x = FALSE, ...) {
      row_means = rowMeans(y) # area of individual ice curves
      ypred = mean(row_means) # area of PDP
      sum((row_means - ypred)^2) # Square difference of area deviation from PDP area for each ICE curve
    }
    input.data = compute_data_for_ice_splitting(effect, testdata = testdata)
  }
  else if (objective == "SS_sd") { # The sample population was divided according to the variance standard deviation predicted by the model
    pdp.feat = effect$features # feature of interest
    split.feats = setdiff(names(testdata), pdp.feat) # features for splitting

    # The ys are the predictions (in this case, the standard deviation)
    X = setDT(testdata)
    Y = setDT(effect$predictor$predict(X))

    split.objective = function(y, x, requires.x = FALSE, ...) {
      y = y$pred
      sum((y - mean(y))^2)
    }
    #split.feats = setdiff(names(testdata), pdp.feat)
    input.data = list(X = X[, ..split.feats, drop = FALSE], Y = Y)  #..split.feats: 引用split.feats变量的值作为列名向量
  }
  else {
    stop(paste("Objective", objective, "is not supported."))
  }

  # Initialize the parent node of the tree
  parent = Node$new(id = 0, depth = 1, subset.idx = seq_len(nrow(input.data$X)), improvement.met = FALSE, intImp = 0)

  # Perform splitting for the parent
  tree = list(list(parent)) # Each layer is a list, and initially the first layer contains the root node.

  for (depth in seq_len(n.split)) { # From layer 1 to the maximum number of layers

    leaves = tree[[depth]] # All nodes in the current layer

    tree[[depth + 1]] = list() # Preparing the container for the next layer

    for (node.idx in seq_along(leaves)) { # Iterate over each node in the current layer

      node.to.split = leaves[[node.idx]]

      if (!is.null(node.to.split)) {
        # Find the best split
        node.to.split$computeSplit(X = input.data$X, Y = input.data$Y,
                                   objective = split.objective, impr.par = impr.par,
                                   optimizer = find_best_binary_split, min.split = min.split)
        # Generate left and right nodes
        node.to.split$computeChildren(input.data$X, input.data$Y, split.objective)

        tree[[depth + 1]] = c(tree[[depth + 1]], node.to.split$children)
      } else {
        tree[[depth + 1]] = c(tree[[depth + 1]], list(NULL,NULL))
      }
    }
  }

  return(tree)
}








#' Prepare ICE Matrix and Split Features for Tree Construction
#'
#' Extracts and centers ICE curves from an \code{iml::FeatureEffect} object, and assembles the feature matrix
#' excluding the feature of interest. Intended for use in building interpretable trees from ICE behavior.
#'
#' @param effect An object of class \code{FeatureEffect} (from the \pkg{iml} package) generated with method = "ice".
#'               It must contain long-format ICE results for a single feature.
#' @param testdata A \code{data.frame} or \code{data.table} used to compute ICE, containing all candidate features.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{X}}{A \code{data.table} of splitting features (excluding the feature of interest).}
#'   \item{\code{Y}}{A centered matrix of ICE values, with each row corresponding to an individual sample.}
#' }
#'
#' @details
#' This function transforms ICE data into a wide-format matrix where each row is an ICE curve,
#' and removes any effects of the global PDP by centering each curve around zero.
#' The splitting features matrix \code{X} excludes the feature used for ICE generation.
#'
#' @seealso \code{\link{compute_tree}}, \code{\link{FeatureEffect}}
#'
#' @export
compute_data_for_ice_splitting = function(effect, testdata) {

  df = setDT(testdata)
  df$.id = seq_row(df) # Unique number .id for each sample for subsequent alignment of ICE curves

  ice.feat = effect$features # ice.feat: feature of interest
  features = names(testdata)

  # Candidate features we consider splitting
  split.feats = setdiff(features, ice.feat)
  df.sub = df[, c(".id", split.feats), with = FALSE] # Retain .id and all candidate feature columns, form split feature matrix

  effectdata = effect$results[[1]] # A long-format data.frame
  effectdata = effectdata[effectdata$.type=="ice",]

  # From long-format → wide-format, now each row is the ICE curve of one sample.
  Y = tidyr::spread(effectdata, .borders, .value)
  Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]

  # center ICE curves by their mean
  Y = Y - rowMeans(Y)
  Y = setDT(Y)


  X = df[, split.feats, with = FALSE]

  return(list(X = X, Y = Y))
}








#' Extract Split Criteria from ICE Tree
#'
#' Extracts all internal (non-terminal) node information from a tree object built using \code{compute_tree()}.
#' This includes each node's split feature, split value, improvement metrics, and objective values.
#' Nodes that are terminal (leaf nodes) are ignored in the returned output, but counted.
#'
#' @param tree A list-of-lists structure representing a binary tree built with class \code{Node}, typically from \code{compute_tree()}.
#'
#' @return A \code{data.frame} containing the following columns for internal nodes:
#' \describe{
#'   \item{depth}{Tree depth of the node}
#'   \item{id}{Node ID (not globally unique)}
#'   \item{objective.value}{Objective value after split}
#'   \item{objective.value.parent}{Objective value before split}
#'   \item{intImp}{Relative improvement compared to parent node}
#'   \item{split.feature}{Feature used for splitting}
#'   \item{split.value}{Split threshold or category}
#'   \item{n.final}{Total number of terminal (leaf) nodes}
#' }
#'
#' @examples
#' \dontrun{
#'   tree <- compute_tree(effect, testdata, n.split = 3)
#'   criteria <- extract_split_criteria(tree)
#'   head(criteria)
#' }
#'
#' @export
extract_split_criteria = function(tree){
  list.split.criteria = lapply(tree, function(depth){
    lapply(depth, function(node){

      if(is.null(node$split.feature)){
        df = data.frame("depth" = "final", "id" = "final",
                        "objective.value" = "final",
                        "objective.value.parent" = "final",
                        "intImp" = "final",
                        "split.feature" = "final",
                        "split.value" = "final")
      }
      else{
        df = data.frame("depth" = node$depth, "id" = node$id,
                        "objective.value" = node$objective.value,
                        "objective.value.parent" = node$objective.value.parent,
                        "intImp" = node$intImp,
                        "split.feature" = node$split.feature,
                        "split.value" = node$split.value)
      }
      df
    })
  })

  list.split.criteria = list.clean(list.split.criteria, function(x) length(x) == 0L, TRUE) # Remove empty list items
  df.split.criteria = unlist(list.split.criteria, recursive = FALSE)
  df.split.criteria = as.data.frame(do.call(rbind, df.split.criteria)) # Flatten the nested list into a data.frame, with each row corresponding to a node.
  n.final = length(which(df.split.criteria$depth == "final"))
  df.split.criteria$n.final = n.final
  df.split.criteria = df.split.criteria[df.split.criteria$depth!="final",] # In the end, only the nodes with real splits (intermediate nodes) are retained for plotting or analysis

  return(df.split.criteria)
}








#------------------------------------------------------------------------------------------------------------
# Functions for plotting






#' Prepare Tree Layout for Plotting
#'
#' Generates a layout data frame containing positions and labels for nodes in a binary tree structure
#' returned by \code{compute_tree()}. Intended for use with downstream plotting functions.
#'
#' @param tree A list-of-lists tree structure (output of \code{compute_tree()}) where each list level corresponds to a tree depth
#'             and contains \code{Node} objects or \code{NULL}.
#'
#' @return A \code{data.frame} with one row per non-null node, including the following columns:
#' \describe{
#'   \item{\code{id}}{String ID combining depth and index (e.g., "2_1").}
#'   \item{\code{node.id}}{Node ID (internal integer identifier).}
#'   \item{\code{depth}}{Integer tree depth level.}
#'   \item{\code{index}}{Node position within level.}
#'   \item{\code{x}, \code{y}}{Coordinates for plotting layout (horizontal and vertical positions).}
#'   \item{\code{label}}{Text label to annotate node, including split feature and threshold, or "Leaf".}
#' }
#'
#' @details
#' This function arranges tree nodes along levels (y-axis) and assigns horizontal positions (x) based on their order in each depth level.
#' Useful as input to \code{plot_tree_structure()} for visualizing the hierarchy.
#'
#' @seealso \code{\link{plot_tree_structure}}, \code{\link{compute_tree}}
#'
#' @export
prepare_tree_layout <- function(tree) {
  layout = data.frame()

  for (depth in seq_along(tree)) {
    nodes = tree[[depth]]
    n.nodes = length(nodes)

    for (i in seq_along(nodes)) {
      node = nodes[[i]]
      if (!is.null(node)) {
        label = if (!is.null(node$split.feature)) {
          paste0(node$split.feature, "\n≤ ", round(node$split.value, 3))
        } else {
          "Leaf"
        }
        layout = rbind(layout, data.frame(
          id = paste0(depth, "_", i),
          node.id = node$id,
          depth = depth,
          index = i,
          x = i,
          y = -depth,
          label = label
        ))
      }
    }
  }
  return(layout)
}




#' Visualize Tree Structure from ICE-Based Model
#'
#' Plots the hierarchical structure of a decision tree built on ICE curves using \code{ggplot2},
#' including node labels and connecting edges. This provides an overview of the tree topology without displaying ICE curves.
#'
#' @param tree A list-of-lists structure returned by \code{compute_tree()}, containing \code{Node} objects organized by depth.
#' @param layout A \code{data.frame} generated by \code{prepare_tree_layout()}, specifying coordinates and labels for each node.
#'
#' @return A \code{ggplot} object showing the tree structure, including:
#' \itemize{
#'   \item Labeled nodes indicating the feature and threshold (or "Leaf").
#'   \item Directed edges linking parent to children.
#'   \item Distinct colors per tree depth.
#' }
#'
#' @details
#' Each node is represented by a labeled box, and arrows indicate parent-child relationships.
#' Colors are automatically assigned to indicate depth level.
#'
#' Useful for debugging or structural analysis of interpretable tree models built using ICE data.
#'
#' @seealso \code{\link{prepare_tree_layout}}, \code{\link{compute_tree}}
#'
#' @export
plot_tree_structure <- function(tree, layout) {
  edges = data.frame()

  for (depth in seq_along(tree)) {
    nodes = tree[[depth]]
    if (depth < length(tree)) {
      for (i in seq_along(nodes)) {
        node = nodes[[i]]
        if (!is.null(node)) {
          left.child = node$children$left.child
          right.child = node$children$right.child

          if (!is.null(left.child)) {
            edges = rbind(edges, data.frame(
              x = layout$x[layout$id == paste0(depth, "_", i)],
              y = layout$y[layout$id == paste0(depth, "_", i)],
              xend = layout$x[layout$id == paste0(depth+1, "_", (i-1)*2+1)],
              yend = layout$y[layout$id == paste0(depth+1, "_", (i-1)*2+1)]
            ))
          }

          if (!is.null(right.child)) {
            edges = rbind(edges, data.frame(
              x = layout$x[layout$id == paste0(depth, "_", i)],
              y = layout$y[layout$id == paste0(depth, "_", i)],
              xend = layout$x[layout$id == paste0(depth+1, "_", (i-1)*2+2)],
              yend = layout$y[layout$id == paste0(depth+1, "_", (i-1)*2+2)]
            ))
          }
        }
      }
    }
  }

  # 给每一层设定颜色
  n.depths = length(unique(layout$depth))
  palette = scales::hue_pal()(n.depths)

  layout$color = palette[as.numeric(as.factor(layout$depth))]

  p = ggplot() +
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.15, "cm")),
                 size = 0.7, color = "grey30") +
    geom_label(data = layout, aes(x = x, y = y, label = label, fill = as.factor(depth)),
               color = "black", label.size = 0.3, label.padding = unit(0.2, "lines")) +
    scale_fill_manual(values = palette) +
    theme_void() +
    theme(legend.position = "none")

  return(p)
}








#' Visualize ICE and PDP for Each Tree Split
#'
#' Creates a faceted plot for each non-terminal split in an ICE-based decision tree,
#' showing both individual ICE curves and the average PDP within each resulting subgroup.
#'
#' @param tree A list-of-lists decision tree as returned by \code{compute_tree()}.
#' @param effect A \code{FeatureEffect} object from \pkg{iml}, containing ICE results for a single feature.
#' @param feature_name Optional character string specifying the feature of interest. If \code{NULL}, inferred from \code{effect}.
#' @param ice.alpha Numeric value between 0 and 1 indicating transparency level of individual ICE lines. Default is 0.3.
#' @param pdp.color Color to use for the PDP lines in each facet. Default is \code{"blue"}.
#' @param return.plots Logical; if \code{TRUE}, returns a named list of \code{ggplot} objects instead of printing them.
#'
#' @return If \code{return.plots = TRUE}, a named list of \code{ggplot} objects, one per parent node split.
#'         Otherwise, the plots are printed and \code{NULL} is returned invisibly.
#'
#' @details
#' Each plot corresponds to one binary split in the ICE-based tree. The parent node’s ICE curves
#' are divided into left and right child groups based on the split feature and threshold.
#' Each facet shows the ICE behavior and PDP of that subgroup.
#'
#' The function uses the ordering in \code{tree} to determine depth and node index. The ICE curves
#' are automatically mean-centered for interpretability.
#'
#' @seealso \code{\link{compute_tree}}, \code{\link{plot_tree_structure}}
#'
#' @export
plot_tree <- function(tree, effect, feature_name = NULL, ice.alpha = 0.3,
                      pdp.color = "blue", return.plots = TRUE) {

  if (is.null(feature_name)) {
    feature_name <- effect$features
    if (length(feature_name) != 1) stop("Please provide a single feature_name.")
  }

  # 提取特征编号用于axis label
  feat_num <- gsub(".*?(\\d+)$", "\\1", feature_name)

  # Prepare ICE data
  effectdata <- effect$results[[1]]
  effectdata <- effectdata[effectdata$.type == "ice", ]
  ice_wide <- spread(effectdata, .borders, .value)
  ice_wide <- ice_wide[, setdiff(colnames(ice_wide), c(".type", ".feature")), drop = FALSE]
  ice_wide <- as.data.table(ice_wide)
  ice_wide <- ice_wide[order(.id)]
  ice_wide_centered <- copy(ice_wide)
  row_means <- rowMeans(ice_wide_centered[, -1, with = FALSE])
  ice_wide_centered[, (2:ncol(ice_wide_centered)) := lapply(.SD, function(x) x - row_means), .SDcols = 2:ncol(ice_wide_centered)]

  # Extract grid points (feature values) from colnames
  feature_grid <- colnames(ice_wide_centered)[-1]
  feature_grid <- as.numeric(feature_grid)

  # Convert to long format
  ice_long <- melt(ice_wide_centered, id.vars = ".id",
                   variable.name = ".feature_grid", value.name = ".value")
  ice_long[, .feature_grid := as.numeric(as.character(.feature_grid))]

  # Align testdata(X) to ICE order
  testdata <- as.data.table(effect$predictor$data$X)
  testdata$.id <- seq_len(nrow(testdata))

  plots <- list()

  for (depth in seq_len(length(tree) - 1)) {
    parents <- tree[[depth]]
    children <- tree[[depth + 1]]

    for (j in seq_along(parents)) {
      parent_node <- parents[[j]]
      left_child <- children[[2 * j - 1]]
      right_child <- children[[2 * j]]

      if (!is.null(parent_node) && (!is.null(left_child) || !is.null(right_child))) {
        plot_data <- list()

        split.feature <- parent_node$split.feature
        split.value <- parent_node$split.value

        if (is.null(split.feature) || is.null(split.value)) next

        # Define split labels for left and right
        left_label <- if (is.numeric(split.value)) {
          paste0(split.feature, " ≤ ", format(split.value, digits = 4))
        } else {
          paste0(split.feature, " = ", split.value)
        }

        right_label <- if (is.numeric(split.value)) {
          paste0(split.feature, " > ", format(split.value, digits = 4))
        } else {
          paste0(split.feature, " != ", split.value)
        }

        n_left <- if (!is.null(left_child)) length(left_child$subset.idx[left_child$subset.idx > 0]) else 0
        n_right <- if (!is.null(right_child)) length(right_child$subset.idx[right_child$subset.idx > 0]) else 0

        if (!is.null(left_child) && n_left > 0) {
          id_sub_left <- left_child$subset.idx[left_child$subset.idx > 0]
          ice_left <- ice_long[.id %in% id_sub_left]
          ice_left$node <- paste0(left_label, "\n(n = ", n_left, ")")
          ice_left$node_order <- 1
        } else {
          ice_left <- NULL
        }

        if (!is.null(right_child) && n_right > 0) {
          id_sub_right <- right_child$subset.idx[right_child$subset.idx > 0]
          ice_right <- ice_long[.id %in% id_sub_right]
          ice_right$node <- paste0(right_label, "\n(n = ", n_right, ")")
          ice_right$node_order <- 2
        } else {
          ice_right <- NULL
        }

        ice_node <- rbindlist(list(ice_left, ice_right), fill = TRUE)

        if (nrow(ice_node) > 0) {
          pdp_node <- ice_node[, .(pdp = mean(.value)), by = .(node, .feature_grid)]

          # Make node an ordered factor to control facet order
          ice_node$node <- factor(ice_node$node, levels = unique(ice_node[order(node_order)]$node))
          pdp_node$node <- factor(pdp_node$node, levels = levels(ice_node$node))

          p <- ggplot(ice_node, aes(x = .feature_grid, y = .value, group = .id)) +
            geom_line(alpha = ice.alpha, color = "grey60") +
            geom_line(data = pdp_node, aes(x = .feature_grid, y = pdp, group = node),
                      inherit.aes = FALSE, color = pdp.color, size = 1.2) +
            facet_wrap(~node) +
            theme_bw(base_size = 14) +
            labs(title = paste0("Split at depth ", depth, ", parent id ", parent_node$id),
                 x = bquote(x[.(feat_num)]),
                 y = bquote(hat(f)[.(feat_num)]^{PD}))

          plots[[paste0("split_depth", depth, "_parent", parent_node$id)]] <- p
        }
      }
    }
  }

  if (return.plots) {
    return(plots)
  } else {
    for (p in plots) print(p)
    invisible(NULL)
  }
}






# Same as `plot_tree()`, with an additional parent node above each pair of nodes after one split
plot_tree_plus <- function(tree, effect, feature_name = NULL, ice.alpha = 0.3,
                      pdp.color = "blue", return.plots = FALSE) {

  if (is.null(feature_name)) {
    feature_name <- effect$features
    if (length(feature_name) != 1) stop("Please provide a single feature_name.")
  }

  feat_num <- gsub(".*?(\\d+)$", "\\1", feature_name)

  effectdata <- effect$results[[1]]
  effectdata <- effectdata[effectdata$.type == "ice", ]
  ice_wide <- spread(effectdata, .borders, .value)
  ice_wide <- ice_wide[, setdiff(colnames(ice_wide), c(".type", ".feature")), drop = FALSE]
  ice_wide <- as.data.table(ice_wide)
  ice_wide <- ice_wide[order(.id)]
  ice_wide_centered <- copy(ice_wide)
  row_means <- rowMeans(ice_wide_centered[, -1, with = FALSE])
  ice_wide_centered[, (2:ncol(ice_wide_centered)) := lapply(.SD, function(x) x - row_means), .SDcols = 2:ncol(ice_wide_centered)]

  feature_grid <- colnames(ice_wide_centered)[-1]
  feature_grid <- as.numeric(feature_grid)

  ice_long <- melt(ice_wide_centered, id.vars = ".id",
                   variable.name = ".feature_grid", value.name = ".value")
  ice_long[, .feature_grid := as.numeric(as.character(.feature_grid))]

  testdata <- as.data.table(effect$predictor$data$X)
  testdata$.id <- seq_len(nrow(testdata))

  plots <- list()

  for (depth in seq_len(length(tree) - 1)) {
    parents <- tree[[depth]]
    children <- tree[[depth + 1]]

    for (j in seq_along(parents)) {
      parent_node <- parents[[j]]
      left_child <- children[[2 * j - 1]]
      right_child <- children[[2 * j]]

      if (!is.null(parent_node) && (!is.null(left_child) || !is.null(right_child))) {
        plot_data <- list()

        split.feature <- parent_node$split.feature
        split.value <- parent_node$split.value

        if (is.null(split.feature) || is.null(split.value)) next

        left_label <- if (is.numeric(split.value)) {
          paste0(split.feature, " ≤ ", format(split.value, digits = 4))
        } else {
          paste0(split.feature, " = ", split.value)
        }

        right_label <- if (is.numeric(split.value)) {
          paste0(split.feature, " > ", format(split.value, digits = 4))
        } else {
          paste0(split.feature, " != ", split.value)
        }

        id_parent <- parent_node$subset.idx[parent_node$subset.idx > 0]
        parent_data <- ice_long[.id %in% id_parent]
        pdp_parent <- parent_data[, .(pdp = mean(.value)), by = .feature_grid]

        parent_plot <- ggplot(parent_data, aes(x = .feature_grid, y = .value, group = .id)) +
          geom_line(alpha = ice.alpha, color = "grey60") +
          geom_line(data = pdp_parent, aes(x = .feature_grid, y = pdp),
                    inherit.aes = FALSE, color = pdp.color, size = 1.2) +
          theme_bw(base_size = 14) +
          theme(strip.background = element_rect(fill = "grey90", color = "black"),
                plot.title = element_text(face = "plain", size = 14, hjust = 0.5)) +
          labs(title = paste0("Parent Node (id = ", parent_node$id, ", n = ", length(id_parent), ")"),
               x = bquote(x[.(feat_num)]),
               y = bquote(hat(f)[.(feat_num)]^{PD}))

        id_sub_left <- if (!is.null(left_child)) left_child$subset.idx[left_child$subset.idx > 0] else integer(0)
        id_sub_right <- if (!is.null(right_child)) right_child$subset.idx[right_child$subset.idx > 0] else integer(0)

        plot_left <- NULL
        if (length(id_sub_left) > 0) {
          data_left <- ice_long[.id %in% id_sub_left]
          pdp_left <- data_left[, .(pdp = mean(.value)), by = .feature_grid]
          plot_left <- ggplot(data_left, aes(x = .feature_grid, y = .value, group = .id)) +
            geom_line(alpha = ice.alpha, color = "grey60") +
            geom_line(data = pdp_left, aes(x = .feature_grid, y = pdp),
                      inherit.aes = FALSE, color = pdp.color, size = 1.2) +
            theme_bw(base_size = 14) +
            theme(strip.background = element_rect(fill = "grey90", color = "black"),
                  plot.title = element_text(face = "plain", size = 14, hjust = 0.5)) +
            labs(title = paste0(left_label, "\n(n = ", length(id_sub_left), ")"),
                 x = bquote(x[.(feat_num)]),
                 y = bquote(hat(f)[.(feat_num)]^{PD}))
        }

        plot_right <- NULL
        if (length(id_sub_right) > 0) {
          data_right <- ice_long[.id %in% id_sub_right]
          pdp_right <- data_right[, .(pdp = mean(.value)), by = .feature_grid]
          plot_right <- ggplot(data_right, aes(x = .feature_grid, y = .value, group = .id)) +
            geom_line(alpha = ice.alpha, color = "grey60") +
            geom_line(data = pdp_right, aes(x = .feature_grid, y = pdp),
                      inherit.aes = FALSE, color = pdp.color, size = 1.2) +
            theme_bw(base_size = 14) +
            theme(strip.background = element_rect(fill = "grey90", color = "black"),
                  plot.title = element_text(face = "plain", size = 14, hjust = 0.5)) +
            labs(title = paste0(right_label, "\n(n = ", length(id_sub_right), ")"),
                 x = bquote(x[.(feat_num)]),
                 y = bquote(hat(f)[.(feat_num)]^{PD}))
        }

        full_plot <- gridExtra::grid.arrange(
          parent_plot,
          gridExtra::arrangeGrob(plot_left, plot_right, ncol = 2),
          ncol = 1,
          heights = c(1, 1.2)
        )

        plots[[paste0("split_depth", depth, "_parent", parent_node$id)]] <- full_plot
      }
    }
  }

  if (!return.plots) {
    for (p in plots) grid::grid.draw(p)
    invisible(NULL)
  } else {
    return(plots)
  }
}
