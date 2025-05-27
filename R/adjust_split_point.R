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
