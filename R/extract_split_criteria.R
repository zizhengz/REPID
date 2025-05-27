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
