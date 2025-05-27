#' Plot Relative Interaction Importance for Each Split
#'
#' Visualizes the relative interaction importance (\code{intImp}) of each split in an ICE-based regression tree.
#' This plot helps identify which splits contribute most to interaction effects in the model.
#'
#' @param tree A tree object constructed by \code{compute_tree()}, containing recursive split information and interaction importance scores.
#' It is expected that the tree includes attributes extractable via \code{extract_split_criteria()}.
#'
#' @return A \code{ggplot2} object showing a bar chart of interaction importance values across tree nodes,
#' with bars colored by tree depth.
#'
#' @details
#' The function uses internal metadata such as split depth, node ID, and feature name to label each bar.
#' Interaction importance values (\code{intImp}) are computed externally and assumed to be attached to each split node.
#'
#' @seealso \code{\link{compute_tree}}, \code{\link{extract_split_criteria}}
#'
#' @export

plot_intImp <- function(tree) {
  criteria <- extract_split_criteria(tree)
  criteria <- as.data.table(criteria)
  criteria[, depth := as.integer(as.character(depth))]
  criteria[, intImp := as.numeric(as.character(intImp))]
  criteria[, node_label := paste0("depth ", depth, ", id ", id, "\n", split.feature)]

  ggplot(criteria, aes(x = reorder(node_label, depth), y = intImp, fill = as.factor(depth))) +
    geom_col(width = 0.7) +
    geom_text(aes(label = round(intImp, 3)), vjust = -0.2, size = 3.5) +
    labs(
      title = "Relative Interaction Importance (intImp) per Split",
      x = "Parent Node (depth, id, split feature)",
      y = "Relative Interaction Importance",
      fill = "depth"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

