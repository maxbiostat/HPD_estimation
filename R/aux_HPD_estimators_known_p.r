get_hpd_intervals_Doss_knownP <- function(X,
                                   ci_level = 0.95,
                                   ps = c(0.025, 0.975)) {
  Z <- qnorm(p = (1 + ci_level) / 2)
  
  BL.ests <- mcmcse::mcse.q(x = X,
                            q = ps[1])
  BU.ests <- mcmcse::mcse.q(x = X,
                            q = ps[2])
  
  L.stuff <- BL.ests$est[1] + Z * c(-1, 1) * BL.ests$se[1]
  U.stuff <- BU.ests$est[1] + Z * c(-1, 1) * BU.ests$se[1]
  
  ans <- tibble::tibble(
    point = c(BL.ests$est, BU.ests$est),
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("HPD_L", "HPD_U"),
    method = "Doss"
  )
  return(ans)
}

get_hpd_intervals_Doss_knownP2 <- function(X,
                                          ci_level = 0.95,
                                          ps = c(0.025, 0.975)) {
  Z <- qnorm(p = (1 + ci_level) / 2)
  
  BL.ests <- mcmcse::mcse.q(x = X,
                            q = ps[1], method = "sub")
  BU.ests <- mcmcse::mcse.q(x = X,
                            q = ps[2], method = "sub")
  
  L.stuff <- BL.ests$est[1] + Z * c(-1, 1) * BL.ests$se[1]
  U.stuff <- BU.ests$est[1] + Z * c(-1, 1) * BU.ests$se[1]
  
  ans <- tibble::tibble(
    point = c(BL.ests$est, BU.ests$est),
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("HPD_L", "HPD_U"),
    method = "Doss_SBM"
  )
  return(ans)
}
##
get_hpd_intervals_Meeker_knownP <- function(X,
                                     ci_level = 0.95,
                                     ps = c(0.025, 0.975)) {
  
  HPD <- quantile(X, probs = ps)
  
  sX <- sort(X)
  L.stuff <- sX[quantile_CI(n = length(X), q = ps[1],
                            alpha = ci_level)$Interval]
  U.stuff <- sX[quantile_CI(n = length(X), q = ps[2],
                            alpha = ci_level)$Interval]
  
  ans <- tibble::tibble(
    point = HPD,
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("HPD_L", "HPD_U"),
    method = "Meeker"
  )
  return(ans)
}

get_hpd_intervals_batches_knownP <- function(X,
                                      b_size = 0.1,
                                      ci_level = 0.95,
                                      ps = c(0.025, 0.975)) {
  B <- round(b_size * length(X))
  batches <- split(X, ceiling(seq_along(X) / B))
  hpds <- do.call(rbind,
                  lapply(batches,
                         quantile, probs = ps))
  qs <- c(1 - ci_level, 1 + ci_level) / 2
  L.int <- quantile(hpds[, 1], probs = qs)
  U.int <- quantile(hpds[, 2], probs = qs)
  ans <- tibble::tibble(
    point = colMeans(hpds),
    lwr = c(L.int[1], U.int[1]),
    upr = c(L.int[2], U.int[2]),
    quantity = c("HPD_L", "HPD_U"),
    method = "Batches"
  )
  return(ans)
}