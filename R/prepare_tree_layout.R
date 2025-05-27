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
