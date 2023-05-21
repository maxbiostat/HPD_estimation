source("aux_fake_MCMC.r")
### Symmetric : (4.01, 0.0352)
### Asymmetric : (2.58, 0.254)
### Highly asymmetric : (0, 1)
## should the target be log-normal (TRUE) or normal (FALSE)?
logNormal <- TRUE
Mu <- 2.58
Sigma <- .254
eff <- .5
Phi <- (1 - eff) / (eff + 1) ## if Phi = 0, independent samples
Alpha <- 0.95
M <- 4E3
X <- generate_fake_MCMC(
  N = M,
  mu = Mu,
  sigma = Sigma,
  phi = Phi,
  LN = logNormal
)

plot(X, type = "l", lwd = 2)
plot(log(X), type = "l", lwd = 2)

coda::effectiveSize(X)
posterior::ess_basic(X)
eff * M
