library(tidyverse)
source("aux.r")
source("aux_fake_MCMC.r")
##############################
############# Functions

run_once <- function(i, Mu, Sigma, Alpha, eff, M, logn) {
  Samples <- generate_fake_MCMC(
    N = M,
    phi = (1 - eff) / (eff + 1),
    mu = Mu,
    sigma = Sigma,
    LN = logn
  )
  #
  HPD <- HDInterval::hdi(Samples, credMass = Alpha)
  BCI <- quantile(Samples, probs = c(1 - Alpha, 1 + Alpha) / 2)
  xbar <- mean(Samples)
  md <- median(Samples)
  #
  qess <- posterior::ess_quantile(x = Samples, probs = c(1 - Alpha, 1 + Alpha) / 2)
  out <- tibble::tibble(
    min = min(Samples),
    max = max(Samples),
    coda_ess = coda::effectiveSize(Samples),
    vats_ess = mcmcse::ess(Samples),
    stan_ess = posterior::ess_basic(Samples),
    tail_ess = posterior::ess_tail(Samples),
    BL_ess = qess[1],
    BU_ess = qess[2],
    point = c(HPD, BCI, xbar, md),
    quantity = c(
      "HPD_L",
      "HPD_U",
      "BCI_L",
      "BCI_U",
      "sample_mean",
      "sample_median"
    ),
    type = "IID",
    target_ESS = round(eff * M),
    M = M,
    replicate = i
  )
  return(out)
}

do_for_ESS <- function(ess) {
  simus <- do.call(rbind,
                   parallel::mclapply(seq_len(Nrep),
                                      function(i) {
                                        run_once(
                                          i,
                                          M = M,
                                          Mu = Mu,
                                          Sigma = Sigma,
                                          Alpha = Alpha,
                                          eff = ess/M,
                                          logn = logNormal
                                        )
                                      }, mc.cores = 6))
  return(simus)
}

### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
## should the target be log-normal (TRUE) or normal (FALSE)
logNormal <- TRUE

Mu <-  0 # 2.58
Sigma <- 1# 0.254

Alpha <- 0.95
M <- 1E4
Nrep <- 500

if (logNormal) {
  true.Mean <- exp(Mu + Sigma ^ 2 / 2)
  true.Median <- exp(Mu)
  true.BCI <- qlnorm(p = c(1 - Alpha, 1 + Alpha) / 2,
                     mean = Mu,
                     sd = Sigma)
  true.HPD <- lognormal_hpd(alpha = Alpha,
                            lmean = Mu,
                            lsd = Sigma)
} else{
  true.Mean <- Mu
  true.Median <- Mu
  true.HPD <- qnorm(p = c(1 - Alpha, 1 + Alpha) / 2,
                    mean = Mu,
                    sd = Sigma)
  true.BCI <- true.HPD
}

true.quants <- tibble::tibble(
  true = c(true.HPD, true.BCI, true.Mean, true.Median),
  quantity = c(
    "HPD_L",
    "HPD_U",
    "BCI_L",
    "BCI_U",
    "sample_mean",
    "sample_median"
  ),
)

ESSs <- c(100, 200, 625, 1000)

results <- do.call(rbind, lapply(ESSs, do_for_ESS))
results.df <- merge(results, true.quants, by = "quantity")

write.csv(results.df, 
          file = paste0("../saved_data/MCMC_lognormal_m=", Mu,
          "_s=", Sigma, "_results.csv"), 
          row.names = FALSE)