# REPID
# Split
# Author: Julia & Giuseppe
# x: [vector] feature to search split
# y: [matrix] n x g matrix containing n ice curves computed on g grid points
# splits: vector of potential split.points to split x, default are unique x values
search_split = function(x, y, splits = unique(x)) { # quantile(x, seq(0, 1, by = 1/100))
  checkmate::assert_vector(x) # 检查x是否为vector，若不是则终止程序
  checkmate::assert_matrix(y) # 检查y是否为matrix，若不是则终止程序
  # x: n x 1 vector with g grid points
  # y: n x g matrix

  # sort feat
  ord = order(x) # output为排好序的vector元素在原vector中的索引
  x = x[ord] # increasing order排序的x_j vector for any j in C
  y = y[ord,] # 按所使用的x_j的大小排好序的ICE curves(already mean-centered)
  N = nrow(y) # number of obs.

  cum_sum_y = Rfast::colCumSums(y) # column-wise cumulative sum for x_j ≤ certain split t
  cum_sum_y2 = Rfast::colCumSums(y^2)

  objective = vapply(splits, function(split) { # for each element of the vector "splits":
    idx = which(x <= split) # Finds indices of x_j that fall into the left partition (x_j ≤ split).
    if (length(idx) == 0 || length(idx) == N) {
      return(Inf) # Invalid split
    }

    N_L = length(idx) # Number of elements of x_j in left partition

    S_L = cum_sum_y[N_L,] # Sum of ice values for each grid point in left partition
    S_R = cum_sum_y[N,] - S_L # Sum of ice values for each grid point in right partition

    SS_L = cum_sum_y2[N_L,] # Sum of squared ice values for each grid point in left partition
    SS_R = cum_sum_y2[N,] - SS_L # Sum of squared ice values for each grid point in right partition

    var_L = sum(SS_L - S_L^2 / N_L)
    var_R = sum(SS_R - S_R^2 / (N - N_L))

    return(var_L + var_R)
  }, FUN.VALUE = NA_real_) # ensures the output is always a numeric vector (NA_real_ is the placeholder for output type)

  best = which.min(objective)
  split.points = splits[best]

  right = min(x[which(x > split.points)])
  left = max(x[which(x <= split.points)])

  data.frame(split.points = (left + right)/2, objective.value = objective[best])
}

best_split = function(X, y) {
  checkmate::assert_data_frame(X)
  res = data.table::rbindlist(lapply(X, function(x) search_split(x, y))) # for each column of X (X must be a list)
  res$feature = colnames(X)
  res$best.split = res$objective.value == min(res$objective.value)
  return(res[, c("feature", "split.points", "objective.value", "best.split")])
}


# Example
library(data.table)
library(ranger)
library(iml)
library(tidyverse)

set.seed(1)

# Simulate Data
n = 5000
x1 = round(runif(n, -1, 1), 1)
x2 = round(runif(n, -1, 1), 3)
x3 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x4 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.7, 0.3))

# noisy vars
x5 = sample(c(0, 1), size = n, replace = TRUE, prob = c(0.5, 0.5))
x6 = rnorm(n, mean = 1, sd = 5)

# target
y = 0.2*x1 - 8*x2 + ifelse(x3 == 0, I(16*x2),0) + ifelse(x1 > 0, I(8*x2),0)
eps = rnorm(n, 0, 0.1*sd(y))
y = y + eps

dat = data.frame(x1, x2, x3, x4, x5, x6, y)
X = dat[, setdiff(colnames(dat), "y")] #only use the columns in dat that are different from "y"

# Fit a random forest model and compute ICE for x2
mod = ranger(y ~ ., data = dat, num.trees = 500)
pred = function(model, newdata) predict(model, newdata)$predictions
model = Predictor$new(mod, data = X, y = dat$y, predict.function = pred)
effect = FeatureEffect$new(model, method = "ice", grid.size = 20, feature = "x2")

# Get mean-centered ICE curves
eff = as.data.table(effect$results)
#subset(eff,.id==1)$.value-mean(subset(eff,.id==1)$.value)
# group by id (5000 obs. -> 5000 id), then take mean of 20 values for each id, and center
eff = as.data.frame(eff[, .value := (.value - mean(.value)), by = c(".type", ".id")])

# Plot ICE curves: WE WANT TO FIND SUBGROUPS SUCH THAT ICE CURVES ARE HOMOGENOUS
ggplot(eff, aes(x = x2, y = .value)) +
  geom_line(aes(group = .id))

# Get ICE values and arrange them in a horizontal matrix
Y = spread(eff, x2, .value)
Y = Y[, setdiff(colnames(Y), c(".type", ".id"))]
Y = as.matrix(Y) # 5000x20

str(X) # contains our feature values
str(Y) # contains ICE values for each grid point

sp = best_split(X, Y)
sp

# Plot regions with ICE curves after best split
split_feat = sp$feature[sp$best.split] # pick the feature to be split
split_val = sp$split.points[sp$best.split] # and pick the split point
breaks = c(min(X[,split_feat]), split_val, max(X[,split_feat]))

split_region = cut(X[,split_feat], breaks = breaks, include.lowest = TRUE)
eff$.split = split_region[effect$results$.id]
# assign the split results to each obs. of X_C (i.e. each ice curve)

ggplot(eff, aes(x = x2, y = .value)) +
  geom_line(aes(group = .id)) + facet_grid(~ .split)

