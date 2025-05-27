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
