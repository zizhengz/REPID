# Load required packages
library(R6)
library(data.table)
library(assertthat)
library(ggplot2)
library(tidyr)
library(scales)

# ---------------------------------------------
# Objective Function (L2 loss)
SS_L2 <- function(y, x, requires.x = FALSE, ...) {
  ypred = colMeans(as.matrix(y))
  sum(t((t(y) - ypred)^2))
}

# ---------------------------------------------
# Find Best Binary Split
# the generic version of search_split() in "Optimal Split using ICE.R"
find_best_binary_split <- function(xval, y, n.splits = 1, min.node.size = 10, objective, ...) {
  assert_choice(n.splits, choices = 1) #  n.splits必须等于1,否则报错 -> binary tree

  # use different split candidates to perform split
  q = generate_split_candidates(xval, n.quantiles = 100, min.node.size = min.node.size)
  splits = vapply(q, FUN = function(i) {
    perform_split(i, xval = xval, y = y, min.node.size = min.node.size,
                  objective = objective, ...)
  }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)
  # select the split point yielding the minimal objective
  best = which.min(splits)

  return(list(split.points = q[best], objective.value = splits[best])) # 返回最佳分割点以及在此取得的最小loss
}

# 1.分位数初步候选
# Generate candidate split points
# xval: Numeric feature vector
# n.quantiles: Number of quantiles for candidate selection
# min.node.size: Minimum number of observations per node
# return A numeric vector of adjusted candidate split points

generate_split_candidates <- function(xval, n.quantiles = NULL, min.node.size = 10) {
  # 检查min.node.size是否是一个整数或“近似整数”,且它的值不能超过(length(xval)-1)/2的下限整值.这可以防止切分点两侧的数据不足
  assert_integerish(min.node.size, lower = 1, upper = floor((length(xval) - 1)/2))

  xval = sort.int(xval) # 把xval排序,为了后续方便找出满足min.node.size要求的中间切分点
  chunk.ind = seq.int(min.node.size + 1, length(xval) - min.node.size, by = min.node.size) # 确保只选中间那段能满足切分条件的点, i.e.左右都有min.node.size个观测值
  xadj = xval[chunk.ind] # 从排序后的xval中提取这些候选切分值

  if (!is.null(n.quantiles)) { # 是否基于分位数生成切分点
    # to speedup we use only quantile values as possible split points
    # qprobs = seq(1/n.quantiles, (n.quantiles - 1)/n.quantiles, by = 1/n.quantiles)
    qprobs = seq(0, 1, by = 1/n.quantiles)
    q = unique(quantile(xadj, qprobs, type = 1)) # "type = 1" -> empirical quantiles; unique保证切分点不会重复
  } else {
    q = unique(xadj)
  }

  q = adjust_split_point(q, xval)  # use a value between two subsequent points
  return(q) # 最终返回的是一组不与与xval中的值重合的经过处理的候选切分点
}

# 2.先去重复+边界调整
# Adjust split points to lie between unique x values
# split.points: Raw candidate split points
# xval: The sorted x values
# return Adjusted split points

adjust_split_point <- function(split.points, xval) {
  q = split.points # split.points: 原始的切分点,可能与xval中的值重合
  x.unique = sort.int(unique(xval)) # 取出唯一值并排序
  ind = which(x.unique %in% q) # 找出与xval唯一值重合的切分点在x.unique中的位置
  ind = ind[ind < length(x.unique)] # 移除最后一个位置的索引(为了后续步骤中的ind+1不能越界)
  # 如果切分点q中某些值不在x.unique中(没被ind找到), 就用min(diff(...))/2作为微调的最小间隔；
  # 否则,对于每个q, 找出它和下一个唯一值之间的距离,取一半
  eps = if (length(ind) != length(q)) min(diff(x.unique))/2 else (x.unique[ind + 1] - x.unique[ind])/2
  q = q + eps # 所有切分点向右微调一点，落在两个观测值之间

  # 修复边界问题:如果微调后q太靠左或太靠右,调整到相邻两个唯一值的中点
  q[q < x.unique[2]] = mean(x.unique[1:2])
  q[q > x.unique[length(x.unique) - 1]] = mean(x.unique[(length(x.unique) - 1):length(x.unique)])
  return(q)
}

# 3.再映射为合法点
# Replace split.points with closest value from xval taking into account min.node.size
# 将给定的一组split.points替换为数据中最接近它们,并且间隔合理的切分点,以尽可能满足min.node.size约束
get_closest_point = function(split.points, xval, min.node.size = 10) {
  xval = sort.int(xval)
  # try to ensure min.node.size between points (is not guaranteed if many duplicated values exist)
  chunk.ind = seq.int(min.node.size + 1, length(xval) - min.node.size, by = min.node.size) # 确保只选中间那段能满足切分条件的点, i.e.左右都有min.node.size个观测值
  xadj = unique(xval[chunk.ind]) # unique(quantile(xval, prob = chunk.ind/length(xval), type = 1))
  # xval = xval[-c(1, length(xval))]
  split.adj = numeric(length(split.points))
  # 遍历所有目标切分点,把它们替换成最接近的xadj值
  for (i in seq_along(split.adj)) {
    d = xadj - split.points[i]
    ind.closest = which.min(abs(d))
    split.adj[i] = xadj[ind.closest]
    xadj = xadj[-ind.closest] # remove already chosen value
  }

  return(sort.int(split.adj))
}

# Performs a single split and measures the objective
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

# Performs a single split and measures the objective (revised version)
perform_split <- function(split.points, xval, y, min.node.size, objective, ..., max.attempts = 5) {
  split.points = sort.int(split.points)
  split.points = get_closest_point(split.points, xval, min.node.size)

  attempt = 0
  while (attempt < max.attempts) {
    attempt <- attempt + 1
    node.number = findInterval(x = xval, split.points, rightmost.closed = TRUE) + 1 #把样本根据切分点划入区间,得到每个样本所属的子节点编号(因为原返回值从0开始编号,所以+1)
    node.size = tabulate(node.number) # tabulate: 统计整数向量x中每个唯一值的频数,返回一个整数向量

    if (min(node.size) >= min.node.size) {
      y.list = split(y, node.number)
      requires.x = formals(objective)[["requires.x"]]
      x.list = if (isTRUE(requires.x)) split(xval, node.number) else NULL

      res = vapply(seq_along(y.list), FUN = function(i) {
        objective(y = y.list[[i]], x = x.list[[i]], ...)
      }, FUN.VALUE = NA_real_, USE.NAMES = FALSE)

      return(sum(res))
    }

    small.nodes = which(node.size < min.node.size) # 找出哪些切分后区间是“太小”的,需要处理
    if (length(small.nodes) == 0 || length(split.points) == 0) break # 如果没有不合理的节点了(small.nodes为空),或者切分点已经删光(split.points为空),就退出while循环
    drop.index = unique(pmax(small.nodes - 1, 1)) # 将所有不合理的节点编号small.nodes映射成它们对应的左侧切分点索引(也就是split.points中的index),并防止0值索引,去重后准备删除
    drop.index = drop.index[drop.index <= length(split.points)] # 要删的split点索引不能超过当前实际的split数量
    split.points = split.points[-drop.index]

    if (length(split.points) == 0) break
  }

  return(Inf)
}


# performs one split
split_parent_node = function(Y, X, n.splits = 1, min.node.size = 10, optimizer,
                             objective, ...) {
  require(data.table)
  assert_data_frame(X)
  assert_integerish(n.splits)
  assert_integerish(min.node.size)
  assert_function(objective, args = c("y", "x", "requires.x"))
  assert_function(optimizer, args = c("xval", "y"))

  # find best split points per feature
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

# Use it after split_parent_node()
generate_node_index = function(Y, X, result) {
  assert_data_table(result)
  # TODO(Done): fix bug if more than one feature have the same best objective
  #feature = unique(result$feature[result$best.split])
  #split.points = unlist(result$split.points[result$best.split])
  idx = which(result$best.split)[1]
  feature = result$feature[idx]
  split.points = unlist(result$split.points[idx])

  if (is.vector(X))
    xval = X else
      xval = X[, feature]

  cuts = c(min(xval), split.points, max(xval))
  sp = cut(xval, breaks = unique(cuts), include.lowest = TRUE)
  #levels(sp) = paste0(feature, " in ", levels(sp))

  return(list(class = sp, index = split(seq_along(xval), sp)))
}

# 'Node' class
# ---------------------------------------------
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

# ---------------------------------------------
# compute single tree based on Class 'Node'
compute_tree = function(effect, testdata, objective = "SS_L2", n.split, impr.par = 0.05, min.split = 10) {

  if (objective == "SS_L1") {
    split.objective = function(y, x, requires.x = FALSE, ...) {
      require(Rfast)
      ypred = colMeans(as.matrix(y)) # 每列求均值 → PDP曲线
      min(t((t(y) - ypred)^2)) # 用一条最接近平均行为的ICE曲线作为代表
    }
    input.data = compute_data_for_ice_splitting(effect, testdata = testdata)
  }
  else if (objective == "SS_L2") {
    split.objective = function(y, x, requires.x = FALSE, ...) {
      ypred = colMeans(as.matrix(y)) # PDP曲线
      sum(t((t(y) - ypred)^2)) # 所有样本的曲线误差平方和
    }
    input.data = compute_data_for_ice_splitting(effect, testdata = testdata)
  }
  else if (objective == "SS_area") {
    split.objective = function(y, x, requires.x = FALSE, ...) {
      row_means = rowMeans(y) # area of individual ice curves
      ypred = mean(row_means) # area of pdp
      sum((row_means - ypred)^2) # 各ICE曲线的面积偏离PDP面积的平方差
    }
    input.data = compute_data_for_ice_splitting(effect, testdata = testdata)
  }
  else if (objective == "SS_sd") { # 按照模型预测的方差standard deviation来划分样本群体
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
  tree = list(list(parent)) # 每一层是一个list,初始是第一层包含根节点

  for (depth in seq_len(n.split)) { # 从第1层到最大层数

    leaves = tree[[depth]] # 当前层的所有节点

    tree[[depth + 1]] = list() # 准备下一层的容器

    for (node.idx in seq_along(leaves)) { # 遍历当前层的每个节点

      node.to.split = leaves[[node.idx]]

      if (!is.null(node.to.split)) {
        # 查找最优划分
        node.to.split$computeSplit(X = input.data$X, Y = input.data$Y,
                                   objective = split.objective, impr.par = impr.par,
                                   optimizer = find_best_binary_split, min.split = min.split)
        # 生成左右子节点
        node.to.split$computeChildren(input.data$X, input.data$Y, split.objective)

        tree[[depth + 1]] = c(tree[[depth + 1]], node.to.split$children)
      } else {
        tree[[depth + 1]] = c(tree[[depth + 1]], list(NULL,NULL))
      }
    }
  }

  return(tree)
}

compute_data_for_ice_splitting = function(effect, testdata) {

  # effect: effect object of IML method FeatureEffect
  # testdata: X
  # Output: A data.frame where each row corresponds to a ice curve

  df = setDT(testdata) # 把testdata转为data.table
  df$.id = seq_row(df) # 为每个样本加上唯一编号.id,用于后续对齐ICE曲线

  ice.feat = effect$features # ice.feat: feature of interest
  features = names(testdata)

  # Features we consider splitting
  split.feats = setdiff(features, ice.feat) # 其余所有用于切分的候选特征
  df.sub = df[, c(".id", split.feats), with = FALSE] # 保留.id 和所有候选特征列,构成split特征矩阵

  effectdata = effect$results[[1]] # 一个长格式的data.frame
  effectdata = effectdata[effectdata$.type=="ice",] # 只保留.type == "ice"的行

  Y = tidyr::spread(effectdata, .borders, .value) # 从长表 → 宽表,每一行是一个样本的ICE曲线
  Y = Y[, setdiff(colnames(Y), c(".type", ".id", ".feature"))]

  # center ICE curves by their mean
  Y = Y - rowMeans(Y)
  Y = setDT(Y)


  X = df[, split.feats, with = FALSE]

  return(list(X = X, Y = Y))
}

extract_split_criteria = function(tree){
  list.split.criteria = lapply(tree, function(depth){ # 遍历每一层
    lapply(depth, function(node){ # 遍历每一层中的每个节点

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

  list.split.criteria = list.clean(list.split.criteria, function(x) length(x) == 0L, TRUE) # 去掉空列表项
  df.split.criteria = unlist(list.split.criteria, recursive = FALSE)
  df.split.criteria = as.data.frame(do.call(rbind, df.split.criteria)) #把嵌套列表展平成data.frame,每一行对应一个node
  n.final = length(which(df.split.criteria$depth == "final"))
  df.split.criteria$n.final = n.final
  df.split.criteria = df.split.criteria[df.split.criteria$depth!="final",] # 最后只保留了真正有分裂的节点(中间节点),用于绘图或分析

  return(df.split.criteria)
}

#------------------------------------------------------------------------------------------------------------
# functions for plotting

# 准备每个节点的位置布局
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


# 绘制树结构
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

plot_tree <- function(tree, effect, feature_name = NULL, ice.alpha = 0.3,
                      pdp.color = "blue", return.plots = FALSE) {
  require(ggplot2)
  require(tidyr)
  require(data.table)

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
