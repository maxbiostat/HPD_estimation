get_hpd_intervals_Doss <- function(X,
                                   h_level = 0.95,
                                   ci_level = 0.95) {
  
  Z <- qnorm(p = (1 + ci_level) / 2)
  
  HPD <- HDInterval::hdi(X, credMass = h_level)
  
  phat.LB <- estimate_p(q = HPD[1], X = X)
  phat.UB <- estimate_p(q = HPD[2], X = X)
  
  BL.ests <- mcmcse::mcse.q(x = X,
                            q = phat.LB)
  BU.ests <- mcmcse::mcse.q(x = X,
                            q = phat.UB)
  
  ## Equation (12) of Doss et al. (2014)
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
get_hpd_intervals_Stan <- function(X,
                                   h_level = 0.95,
                                   ci_level = 0.95) {
  
  Z <- qnorm(p = (1 + ci_level) / 2)
  
  HPD <- HDInterval::hdi(X, credMass = h_level)
  
  phat.LB <- estimate_p(q = HPD[1], X = X)
  phat.UB <- estimate_p(q = HPD[2], X = X)
  
  SEs <- posterior::mcse_quantile(X, probs = c(phat.LB, phat.UB))
  
  ## Equation (12) of Doss et al. (2014)
  L.stuff <- HPD[1] + Z * c(-1, 1) * SEs[1]
  U.stuff <- HPD[2] + Z * c(-1, 1) * SEs[2]
  
  ans <- tibble::tibble(
    point = HPD,
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("HPD_L", "HPD_U"),
    method = "Stan"
  )
  return(ans)
}
get_hpd_intervals_Doss_SBM <- function(X,
                                   h_level = 0.95,
                                   ci_level = 0.95) {
  Z <- qnorm(p = (1 + ci_level) / 2)
  
  HPD <- HDInterval::hdi(X, credMass = h_level)
  
  phat.LB <- estimate_p(q = HPD[1], X = X)
  phat.UB <- estimate_p(q = HPD[2], X = X)
  
  BL.ests <- mcmcse::mcse.q(x = X,
                            q = phat.LB, method = "sub")
  BU.ests <- mcmcse::mcse.q(x = X,
                            q = phat.UB, method = "sub")
  
  ## Equation (12) of Doss et al. (2014)
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
get_hpd_intervals_Meeker <- function(X,
                                     h_level = 0.95,
                                     ci_level = 0.95) {
  
  HPD <- HDInterval::hdi(X, credMass = h_level)
  
  phat.LB <- estimate_p(q = HPD[1], X = X)
  phat.UB <- estimate_p(q = HPD[2], X = X)
  
  sX <- sort(X)
  L.stuff <- sX[quantile_CI(n = length(X), q = phat.LB,
                            alpha = ci_level)$Interval]
  U.stuff <- sX[quantile_CI(n = length(X), q = phat.UB,
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
get_hpd_intervals_SPIn <- function(X,
                                   h_level = 0.95,
                                   ci_level = 0.95, lb  = 0) {
  
  HPD <- SPIn::SPIn(X, conf = h_level, lb = lb)$spin
  
  phat.LB <- estimate_p(q = HPD[1], X = X)
  phat.UB <- estimate_p(q = HPD[2], X = X)
  
  sX <- sort(X)
  L.stuff <- sX[quantile_CI(n = length(X), q = phat.LB,
                            alpha = ci_level)$Interval]
  U.stuff <- sX[quantile_CI(n = length(X), q = phat.UB,
                            alpha = ci_level)$Interval]
  
  ans <- tibble::tibble(
    point = HPD,
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("HPD_L", "HPD_U"),
    method = "SPIn"
  )
  return(ans)
}
get_hpd_intervals_batches <- function(X,
                                      b_size = 0.1,
                                      h_level = 0.95,
                                      ci_level = 0.95) {
  B <- round(b_size * length(X))
  batches <- split(X, ceiling(seq_along(X) / B))
  hpds <- do.call(rbind,
                  lapply(batches,
                         HDInterval::hdi, credMass = h_level))
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

##
opt_hpd <- function(samples, alpha) {
  opt_int <- function(q1, a) {
    q2 <- q1 + a
    Lest <- mcmcse::mcse.q(x = samples, q = q1)
    Uest <- mcmcse::mcse.q(x = samples, q = q2)
    cand <- c(Lest$est, Uest$est)
    return(cand[2] - cand[1])
  }
  Opt <- optimise(opt_int,
                  a = alpha,
                  lower = 1E-10,
                  upper = 1 - alpha)
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  return(list(
    pL = q1.opt,
    pU = q2.opt,
    Lhat = mcmcse::mcse.q(x = samples, q = q1.opt),
    Uhat = mcmcse::mcse.q(x = samples, q = q2.opt)
  ))
}
###
get_hpd_intervals_DossOpt <- function(X,
                                      h_level = 0.95,
                                      ci_level = 0.95) {
  Z <- qnorm(p = (1 + ci_level) / 2)
  
  OptStuff <- opt_hpd(samples = X, a = h_level)
  
  L.stuff <- OptStuff$Lhat$est + Z * c(-1, 1) * OptStuff$Lhat$se
  U.stuff <- OptStuff$Uhat$est + Z * c(-1, 1) * OptStuff$Uhat$se
  
  ans <- tibble::tibble(
    point = c(OptStuff$Lhat$est, OptStuff$Uhat$est),
    lwr = c(ifelse(is.na(L.stuff[1]), -Inf, L.stuff[1]),
            U.stuff[1]),
    upr = c(L.stuff[2],
            ifelse(is.na(U.stuff[2]), Inf, U.stuff[2])),
    quantity = c("HPD_L", "HPD_U"),
    method = "Doss_Opt"
  )
  return(ans)
}