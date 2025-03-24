parse_bootCI <- function(bt.ci, varname = "") {
  
  normal <- tibble::tibble(
    point = bt.ci$t0,
    lwr = bt.ci$normal[2],
    upr = bt.ci$normal[3],
    phat_lwr = NA,
    phat_upr = NA,
    quantity = varname,
    type = "normal",
    method =  "bootstrap"
  )
  
  basic <- tibble::tibble(
    point = bt.ci$t0,
    lwr = bt.ci$basic[4],
    upr = bt.ci$basic[5],
    phat_lwr = NA,
    phat_upr = NA,
    quantity = varname,
    type = "basic",
    method =  "bootstrap"
  )
  
  percentile <- tibble::tibble(
    point = bt.ci$t0,
    lwr = bt.ci$percent[4],
    upr = bt.ci$percent[5],
    phat_lwr = NA,
    phat_upr = NA,
    quantity = varname,
    type = "percent",
    method =  "bootstrap"
  )
  
  ans <- do.call(rbind,
                 list(normal, percentile, basic))
  
  return(ans)
}

hpd_intervals_boot <- function(samples,
                               alpha = .95,
                               gamma = .95,
                               B = 500,
                               nu = NULL,
                               cores = 10) {
  
  M <- length(samples)
  
  if (is.null(nu)) {
    bsize <- mcmcse::batchSize(x = samples)
    nu <- "estimated"
  } else{
    bsize <- round(M ^ nu)
  }
  
  if (bsize <= 5) {
    HPD.fun.1 <- function(tsb, ind) {
      c(HDInterval::hdi(tsb[ind], credMass = alpha))
    }
    bb <- boot::boot(samples, HPD.fun.1, R = B,
                     parallel = "multicore", ncpus = cores)
  } else{
    HPD.fun.2 <- function(tsb) {
      c(HDInterval::hdi(tsb, credMass = alpha))
    }
    bb <- boot::tsboot(samples,
                       HPD.fun.2,
                       R = B,
                       l = bsize,
                       sim = "fixed",
                       parallel = "multicore", ncpus = cores)
  }
  
  raw.ci.a <- boot::boot.ci(boot.out = bb, conf = gamma,
                            type = c("norm",
                                     "perc", "basic"),
                          index = 1)
  
  raw.ci.b <- boot::boot.ci(boot.out = bb, conf = gamma,
                            type = c("norm",
                                     "perc", "basic"),
                            index = 2)
  
    
  ans.a <- parse_bootCI(raw.ci.a, varname = "HPD_L")
  ans.b <- parse_bootCI(raw.ci.b, varname = "HPD_U")
  
  ans <- rbind(ans.a, ans.b)
  ans$n_boot <- B
  ans$batch_size <- bsize
  ans$nu <- as.character(nu)
  
  return(ans)
}