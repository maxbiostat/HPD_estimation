source("aux_fake_MCMC.r")
source("aux.r")
library(mcmcse)
### Symmetric : (4.01, 0.0352)
### Asymmetric : (2.58, 0.254)
### Highly asymmetric : (0, 1)
## should the target be log-normal (TRUE) or normal (FALSE)?
logNormal <- TRUE
Mu <- 0 #2.58
Sigma <- 1 # .254
eff <- .999
Phi <- (1 - eff) / (eff + 1) ## if Phi = 0, independent samples
Alpha <- 0.95
M <- 1E4
X <- generate_fake_MCMC(
  N = M,
  mu = Mu,
  sigma = Sigma,
  phi = Phi,
  LN = logNormal
)

coda::effectiveSize(X)
posterior::ess_basic(X)
posterior::ess_quantile(X)
eff * M
####
Alpha <- 0.95
true.HPD <- lognormal_hpd(alpha = Alpha, lmean = Mu, lsd = Sigma)
true.ps <- lognormal_hpd_percentiles(alpha = Alpha, lmean = Mu, lsd = Sigma)

range(X)
mcmcse::mcse.q(x = X, q = true.ps[1],
               warn = TRUE)
mcmcse::mcse.q(x = X, q = true.ps[1],
               method = "obm",
               warn = TRUE)
mcmcse::mcse.q(x = X, q = true.ps[2],
               warn = TRUE)
mcmcse::mcse.q(x = X, q = true.ps[2],
               method = "obm",
               warn = TRUE)
est.HPD <- HDInterval::hdi(X)
est.ps <- sapply(est.HPD, function(x) mean(X <= x))

cbind(true.ps, est.ps)
cbind(true.HPD, est.HPD)

source("aux_quantile_CI.r")
source("aux_HPD_estimators.r")
source("aux_HPD_estimators_known_p.r")
get_hpd_intervals_Meeker(X)
get_hpd_intervals_Meeker_knownP(X, ps = true.ps )
get_hpd_intervals_Doss(X)
get_hpd_intervals_DossOpt(X)
get_hpd_intervals_Doss_knownP(X, ps = true.ps)
true.HPD
# get_hpd_intervals_SPIn(X)

mcmcse::mcse.q(X, q = true.ps[1], method = "bm")
mcmcse::mcse.q(X, q = true.ps[1], method = "obm")
mcmcse::mcse.q(X, q = true.ps[1], method = "sub")

mcmcse::mcse.q(X, q = true.ps[2], method = "bm")
mcmcse::mcse.q(X, q = true.ps[2], method = "obm")
mcmcse::mcse.q(X, q = true.ps[2], method = "sub")
