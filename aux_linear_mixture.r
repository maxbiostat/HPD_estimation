parse_fit <- function(fit, type = c("gamma", "norm")) {
  if (type == "norm") {
    res <- list(
      mix_weights = as.numeric(fit[1,]),
      mix_means = as.numeric(fit[2,]),
      mix_sds = as.numeric(fit[3,])
    )
  } else{
    res <- list(
      mix_weights = as.numeric(fit[1,]),
      mix_shapes = as.numeric(fit[2,]),
      mix_rates = as.numeric(fit[3,])
    )
  }
  return(res)
}

mix_dens_gaussian <- function(x, mus, sigmas, ws, log = FALSE) {
  K <- length(mus)
  if (length(sigmas) != K)
    stop("vectors of parameters not the same size")
  if (length(ws) != K)
    stop("vector of weights is not the correct size")
  lws <- log(ws)
  compute_gauss_dens <- function(x) {
    require(matrixStats)
    lps <- rep(NA, K)
    for (k in 1:K) {
      lps[k] <- lws[k] + dnorm(
        x = x,
        mean = mus[k],
        sd = sigmas[k],
        log = TRUE
      )
    }
    return(matrixStats::logSumExp(lps))
  }
  ans <- sapply(x, compute_gauss_dens)
  if (!log)
    ans <- exp(ans)
  return(ans)
}

mix_dens_gamma <- function(x, as, bs, ws, log = FALSE) {
  K <- length(as)
  if (length(bs) != K)
    stop("vectors of parameters not the same size")
  if (length(ws) != K)
    stop("vector of weights is not the correct size")
  lws <- log(ws)
  compute_gamma_dens <- function(x) {
    require(matrixStats)
    lps <- rep(NA, K)
    for (k in 1:K)
      lps[k] <- lws[k] + dgamma(
        x = x,
        shape = as[k],
        rate = bs[k],
        log = TRUE
      )
    return(matrixStats::logSumExp(lps))
  }
  ans <- sapply(x, compute_gamma_dens)
  if (!log)
    ans <- exp(ans)
  return(ans)
}

mix_cdf_gaussian <- function(q, mus, sigmas, ws, log = FALSE) {
  K <- length(mus)
  if (length(sigmas) != K)
    stop("vectors of parameters not the same size")
  if (length(ws) != K)
    stop("vector of weights is not the correct size")
  lws <- log(ws)
  compute_gauss_lcdf <- function(x) {
    require(matrixStats)
    lps <- rep(NA, K)
    for (k in 1:K)
      lps[k] <- lws[k] + pnorm(
        q = x,
        mean = mus[k],
        sd = sigmas[k],
        log.p =  TRUE
      )
    return(matrixStats::logSumExp(lps))
  }
  ans <- sapply(q, compute_gauss_lcdf)
  if (!log)
    ans <- exp(ans)
  return(ans)
}

mix_cdf_gamma <- function(q, as, bs, ws, log = FALSE) {
  K <- length(as)
  if (length(bs) != K)
    stop("vectors of parameters not the same size")
  if (length(ws) != K)
    stop("vector of weights is not the correct size")
  lws <- log(ws)
  compute_gamma_lcdf <- function(x) {
    require(matrixStats)
    lps <- rep(NA, K)
    for (k in 1:K)
      lps[k] <- lws[k] + pgamma(
        q = x,
        shape = as[k],
        rate = bs[k],
        log.p = TRUE
      )
    return(matrixStats::logSumExp(lps))
  }
  ans <- sapply(q, compute_gamma_lcdf)
  if (!log)
    ans <- exp(ans)
  return(ans)
}

mix_invcdf_gaussian <- function(p, mus, sigmas, ws) {
  obj.fn <- function(q) {
    lF <- mix_cdf_gaussian(
      q = q,
      mus = mus,
      sigmas = sigmas,
      ws = ws,
      log = TRUE
    )
    return((lF - log(p)) ^ 2)
  }
  obj.fn <- Vectorize(obj.fn)
  Opt <- optimise(obj.fn,
                  interval = c(min(mus) - 5 * max(sigmas),
                               max(mus) + 5 * max(sigmas)),
                  tol = 1E-20)
  return(as.numeric(Opt$minimum))
}

mix_invcdf_gamma <- function(p, as, bs, ws) {
  obj.fn <- function(q) {
    lF <- mix_cdf_gamma(
      q = q,
      as = as,
      bs = bs,
      ws = ws,
      log = TRUE
    )
    return((lF - log(p)) ^ 2)
  }
  obj.fn <- Vectorize(obj.fn)
  
  Opt <- optimise(obj.fn,
                  interval = c(0,
                               # mean + 5*sd
                               max(as / bs) + 5 * max(sqrt(as) / bs)),
                  tol = 1E-20)
  return(as.numeric(Opt$minimum))
}

mix_hpd_gaussian <- function(alpha, mus, sigmas, ws) {
  opt_int <- function(q1, a) {
    q2 <- q1 + a
    cand <- sapply(c(q1, q2),
                   function(p)
                     mix_invcdf_gaussian(
                       p = p,
                       mus = mus,
                       sigmas = sigmas,
                       ws = ws
                     ))
    return(cand[2] - cand[1])
  }
  Opt <- optimise(opt_int,
                  a = alpha,
                  lower = 1E-10,
                  upper = 1 - alpha)
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  ans <- sapply(c(q1.opt, q2.opt),
                function(p)
                  mix_invcdf_gaussian(
                    p = p,
                    mus = mus,
                    sigmas = sigmas,
                    ws = ws
                  ))
  return(ans)
}

mix_hpd_gamma <- function(alpha, as, bs, ws) {
  opt_int <- function(q1, a) {
    q2 <- q1 + a
    cand <- sapply(c(q1, q2),
                   function(p)
                     mix_invcdf_gamma(
                       p = p,
                       as = as,
                       bs = bs,
                       ws = ws
                     ))
    return(cand[2] - cand[1])
  }
  Opt <- optimise(opt_int,
                  a = alpha,
                  lower = 1E-10,
                  upper = 1 - alpha)
  q1.opt <- Opt$minimum
  q2.opt <- q1.opt + alpha
  ans <- sapply(c(q1.opt, q2.opt),
                function(p)
                  mix_invcdf_gamma(
                    p = p,
                    as = as,
                    bs = bs,
                    ws = ws
                  ))
  return(ans)
}

HPD_CI_gaussian <- function(X,
                            alpha = 0.95,
                            ci.level = 0.95) {
  mfit <- RBesT::automixfit(X)
  pars <- parse_fit(mfit, type = "norm")
  mixf <- function(x)
    mix_dens_gaussian(
      x = x,
      mus = pars$mix_means,
      sigmas = pars$mix_sds,
      ws = pars$mix_weights
    )
  mixf <- Vectorize(mixf)
  HPD.hat <- mix_hpd_gaussian(
    alpha = alpha,
    mus = pars$mix_means,
    sigmas = pars$mix_sds,
    ws = pars$mix_weights
  )
  hat.ps <- mix_cdf_gaussian(
    q = HPD.hat,
    mus = pars$mix_means,
    sigmas = pars$mix_sds,
    ws = pars$mix_weights
  )
  S <- length(X)
  Z <- qnorm(p = (1 + ci.level) / 2)
  se.L <-
    sqrt(exp(log(hat.ps[1]) + log1p(-hat.ps[1])) / (S * (mixf(HPD.hat[1])) ^
                                                      2))
  L.ci <- HPD.hat[1] + c(-1, 1) * Z * se.L
  se.U <-
    sqrt(exp(log(hat.ps[2]) + log1p(-hat.ps[2])) / (S * (mixf(HPD.hat[2])) ^
                                                      2))
  U.ci <- HPD.hat[2] + c(-1, 1) * Z * se.U
  ans <- tibble::tibble(
    point = HPD.hat,
    lwr = c(L.ci[1], U.ci[1]),
    upr = c(L.ci[2], U.ci[2]),
    quantity = c("HPD_L", "HPD_U"),
    method = "Gaussian_mixture"
  )
  return(ans)
}


HPD_CI_gamma <- function(X,
                         alpha = 0.95,
                         ci.level = 0.95) {
  mfit <- RBesT::automixfit(X, type = "gamma")
  pars <- parse_fit(mfit, type = "gamma")
  mixf <- function(x)
    mix_dens_gamma(
      x = x,
      as = pars$mix_shapes,
      bs = pars$mix_rates,
      ws = pars$mix_weights
    )
  mixf <- Vectorize(mixf)
  HPD.hat <- mix_hpd_gamma(
    alpha = alpha,
    as = pars$mix_shapes,
    bs = pars$mix_rates,
    ws = pars$mix_weights
  )
  hat.ps <- mix_cdf_gamma(
    q = HPD.hat,
    as = pars$mix_shapes,
    bs = pars$mix_rates,
    ws = pars$mix_weights
  )
  S <- length(X)
  Z <- qnorm(p = (1 + ci.level) / 2)
  se.L <-
    sqrt(exp(log(hat.ps[1]) + log1p(-hat.ps[1])) /
           (S * (mixf(HPD.hat[1])) ^ 2))
  L.ci <- HPD.hat[1] + c(-1, 1) * Z * se.L
  se.U <- sqrt(exp(log(hat.ps[2]) + log1p(-hat.ps[2])) /
                 (S * (mixf(HPD.hat[2])) ^ 2))
  U.ci <- HPD.hat[2] + c(-1, 1) * Z * se.U
  ans <- tibble::tibble(
    point = HPD.hat,
    lwr = c(L.ci[1], U.ci[1]),
    upr = c(L.ci[2], U.ci[2]),
    quantity = c("HPD_L", "HPD_U"),
    method = "Gamma_mixture"
  )
  return(ans)
}
