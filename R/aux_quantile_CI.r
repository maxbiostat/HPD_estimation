quantile_CI <- function(n, q, alpha = 0.95) {
  # Copied from https://stats.stackexchange.com/questions/99829/how-to-obtain-a-confidence-interval-for-a-percentile
  # Search over a small range of upper and lower order statistics for the
  # closest coverage to 1-alpha (but not less than it, if possible).
  gg <- 1 - alpha
  u <- qbinom(1 - gg / 2, n, q) + (-2:2) + 1
  l <- qbinom(gg / 2, n, q) + (-2:2)
  u[u > n] <- Inf
  l[l < 0] <- -Inf
  coverage <-
    outer(l, u, function(a, b)
      pbinom(b - 1, n, q) - pbinom(a - 1, n, q))
  if (max(coverage) < 1 - gg)
    i <- which(coverage == max(coverage))
  else
    i <- which(coverage == min(coverage[coverage >= 1 - gg]))
  i <- i[1]
  # Return the order statistics and the actual coverage.
  u <- rep(u, each = 5)[i]
  l <- rep(l, 5)[i]
  return(list(Interval = c(l, u),
              Coverage = coverage[i]))
}