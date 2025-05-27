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
