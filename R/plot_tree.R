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
          labs(title = paste0("Parent Node (depth = ", parent_node$depth, ", id = ", parent_node$id, ", n = ", length(id_parent), ")"),
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
