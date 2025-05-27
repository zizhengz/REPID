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
