#' Computes confidence intervals for the BCI endpoints using the method of Doss et al. 2014
#' 
#'
#' @param X vector containing MCMC draws
#' @param b_level probability level for the BCI
#' @param ci_level 
#'
#' @return
#' @export
#'
#' @examples
get_bci_intervals_Doss <- function(X,
                                   b_level = 0.95,
                                   ci_level = 0.95) {
  Z <- qnorm(p = (1 + ci_level) / 2)
  
  BL.ests <- mcmcse::mcse.q(x = X,
                            q = (1 - b_level) / 2)
  BU.ests <- mcmcse::mcse.q(x = X,
                            q = (1 + b_level) / 2)
  
  L.stuff <- BL.ests$est[1] + Z * c(-1, 1) * BL.ests$se[1]
  U.stuff <- BU.ests$est[1] + Z * c(-1, 1) * BU.ests$se[1]
  
  ans <- tibble::tibble(
    point = c(BL.ests$est, BU.ests$est),
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("BCI_L", "BCI_U"),
    method = "Doss"
  )
  return(ans)
}
##
compute_bci  <- function(samples, level = 0.95) {
  stats::quantile(samples, probs = c(1 - level, 1 + level) / 2)
}
##
get_bci_intervals_Meeker <- function(X,
                                     b_level = 0.95,
                                     ci_level = 0.95) {
  BCI <- compute_bci(samples = X, level = b_level)
  
  sX <- sort(X)
  L.stuff <- sX[quantile_CI(n = length(X),
                            q = (1 - b_level) / 2,
                            alpha = ci_level)$Interval]
  U.stuff <- sX[quantile_CI(n = length(X),
                            q = (1 + b_level) / 2,
                            alpha = ci_level)$Interval]
  
  ans <- tibble::tibble(
    point = BCI,
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("BCI_L", "BCI_U"),
    method = "Meeker"
  )
  return(ans)
}
get_bci_intervals_Stan <- function(X,
                                   b_level = 0.95,
                                   ci_level = 0.95) {
  BCI <- compute_bci(samples = X, level = b_level)
  
  Z <- qnorm(p = (1 + ci_level) / 2)
  qs <- (1 + c(-1, 1) * b_level) / 2
  sEs <- posterior::mcse_quantile(X, probs = qs)
  
  L.stuff <- BCI[1] + Z * c(-1, 1) * sEs[1]
  U.stuff <- BCI[2] + Z * c(-1, 1) * sEs[2]
  
  ans <- tibble::tibble(
    point = BCI,
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("BCI_L", "BCI_U"),
    method = "Stan"
  )
  return(ans)
}

get_CLT_se <- function(x, p, f, n) {
  sqrt(p * (1 - p) / (n * f(x) ^ 2))
}
get_bci_intervals_CLT <- function(X,
                                  b_level = 0.95,
                                  ci_level = 0.95) {
  M <- length(X)
  BCI <- compute_bci(samples = X, level = b_level)
  
  Z <- qnorm(p = (1 + ci_level) / 2)
  
  qs <- (1 + c(-1, 1) * b_level) / 2
  dd <- density(X)
  fhat <- approxfun(x = dd$x, y = dd$y)
  sEs <- c(get_CLT_se(
    x = BCI[1],
    p = qs[1],
    f = fhat,
    n = M
  ),
  get_CLT_se(
    x = BCI[2],
    p = qs[2],
    f = fhat,
    n = M
  ))
  
  L.stuff <- BCI[1] + Z * c(-1, 1) * sEs[1]
  U.stuff <- BCI[2] + Z * c(-1, 1) * sEs[2]
  
  ans <- tibble::tibble(
    point = BCI,
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("BCI_L", "BCI_U"),
    method = "CLT"
  )
  return(ans)
}

get_bci_intervals_batches <- function(X,
                                      b_size = 0.1,
                                      b_level = 0.95,
                                      ci_level = 0.95) {
  B <- round(b_size * length(X))
  batches <- split(X, ceiling(seq_along(X) / B))
  bcis <- do.call(rbind,
                  lapply(batches, compute_bci, level = b_level))
  qs <- c(1 - ci_level, 1 + ci_level) / 2
  L.int <- quantile(bcis[, 1], probs = qs)
  U.int <- quantile(bcis[, 2], probs = qs)
  ans <- tibble::tibble(
    point = colMeans(bcis),
    lwr = c(L.int[1], U.int[1]),
    upr = c(L.int[2], U.int[2]),
    quantity = c("BCI_L", "BCI_U"),
    method = "Batches"
  )
  return(ans)
}