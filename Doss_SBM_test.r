library(tidyverse)
source("../aux.r")
source("../aux_fake_MCMC.r")
source("../aux_linear_mixture.r")
source("../aux_estimation_quality.r")
source("../aux_HPD_estimators.r")
##############################
############# Functions
run_once <- function(i, Mu, Sigma, Phi, Alpha, Omega, M, logn) {
  
  Samples <- generate_fake_MCMC(
    N = M,
    mu = Mu,
    sigma = Sigma,
    phi = Phi,
    LN = logn
  )
  
  Doss.ests <- get_hpd_intervals_Doss(X = Samples,
                                      h_level = Alpha,
                                      ci_level = Omega)
  
  
  Doss.sbm.ests <- get_hpd_intervals_Doss_SBM(X = Samples,
                                              h_level = Alpha,
                                              ci_level = Omega)
  Opt.ests <-
    get_hpd_intervals_DossOpt(X = Samples,
                              h_level = Alpha,
                              ci_level = Omega)
  
  concat <- do.call(
    rbind,
    list(
      Doss.ests,
      Doss.sbm.ests,
      Opt.ests
    )
  )
  
  out <- tibble::tibble(
    min = min(Samples),
    max = max(Samples),
    concat,
    replicate = i
  )
  return(out)
}
##########################################
### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
## should the target be log-normal (TRUE) or normal (FALSE)?
logNormal <- TRUE

Mu <- 2.58
Sigma <- .254
eff <- .999
Phi <- (1 - eff) / (eff + 1) ## if Phi = 0, independent samples

Alpha <- 0.95
TargetCoverage <- 0.95
M <- 1E4
Nrep <- 500

if (logNormal) {
  true.HPD <- lognormal_hpd(alpha = Alpha,
                            lmean = Mu,
                            lsd = Sigma)
} else{
  true.HPD <- qnorm(p = c(1 - Alpha, 1 + Alpha) / 2,
                    mean = Mu,
                    sd = Sigma)
}

Ncores <- 10
simu.time <- system.time(
  simus <- do.call(rbind,
                   parallel::mclapply(seq_len(Nrep),
                                      function(i) {
                                        run_once(
                                          i,
                                          M = M,
                                          Mu = Mu,
                                          Sigma = Sigma,
                                          Phi = Phi,
                                          Alpha = Alpha,
                                          Omega = TargetCoverage,
                                          logn = logNormal
                                        )
                                      }, mc.cores = Ncores))
)
simu.time

true.vals <- tibble::tibble(true_value = true.HPD,
                            quantity = c("HPD_L", "HPD_U"))
simus.df <- merge(simus, true.vals, by = "quantity")
simus.df <-
  simus.df %>% mutate(covers = is_in(true_value, lwr, upr))
simus.df <-
  simus.df %>% mutate(lwr = ifelse(is.infinite(lwr), 0, lwr))

coverage <-
  aggregate(covers ~ quantity + method, simus.df, mean) ## coverage
ci.width <-
  aggregate((upr - lwr) ~ quantity + method, simus.df, mean) ## average CI width
bias <-
  aggregate((point - true_value) ~ quantity + method, simus.df, mean)
variance <- aggregate(point ~ quantity + method, simus.df, var)
mse <-
  aggregate((point - true_value) ^ 2 ~ quantity + method, simus.df, mean)
ans <-
  Reduce(function(...)
    merge(..., by = c("quantity", "method")),
    list(coverage, ci.width,
         bias,
         variance, mse))
names(ans) <- c("quantity",
                "method",
                "est_coverage",
                "ci_width",
                "bias",
                "variance",
                "MSE")

ans

ggplot() +
  geom_pointrange(
    data = subset(simus.df, quantity == "HPD_L"),
    mapping = aes(
      x = replicate,
      y = point,
      ymin = lwr,
      ymax = upr,
      colour = covers
    )
  ) +
  geom_hline(yintercept = true.HPD[1], linetype = "longdash") +
  facet_grid(method~., scales = "free_y") +
  ggtitle("HPD Lower") + 
  theme_bw(base_size = 20)

ggplot() +
  geom_pointrange(
    data = subset(simus.df, quantity == "HPD_U"),
    mapping = aes(
      x = replicate,
      y = point,
      ymin = lwr,
      ymax = upr,
      colour = covers
    )
  ) +
  geom_hline(yintercept = true.HPD[2], linetype = "longdash") +
  facet_grid(method~., scales = "free_y") +
  ggtitle("HPD Upper") + 
  theme_bw(base_size = 20)
