library(tidyverse)
source("aux.r")
source("aux_fake_MCMC.r")
source("aux_linear_mixture.r")
source("aux_quantile_CI.r")
source("aux_HPD_estimators.r")
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
  
  Meeker.ests <- get_hpd_intervals_Meeker(X = Samples,
                                          h_level = Alpha,
                                          ci_level = Omega)
  
  SPIn.ests <- get_hpd_intervals_SPIn(X = Samples,
                                      h_level = Alpha,
                                      ci_level = Omega)
  
  Doss.ests <- get_hpd_intervals_Doss(X = Samples,
                                      h_level = Alpha,
                                      ci_level = Omega)
  
  Stan.ests <- get_hpd_intervals_Stan(X = Samples,
                                      h_level = Alpha,
                                      ci_level = Omega)
  
  Doss.sbm.ests <- get_hpd_intervals_Doss_SBM(X = Samples,
                                              h_level = Alpha,
                                              ci_level = Omega)
  
  Mixture.N.ests <- HPD_CI_gaussian(X = Samples,
                                    alpha = Alpha,
                                    ci.level = Omega)
  Mixture.G.ests <- HPD_CI_gamma(X = Samples,
                                 alpha = Alpha,
                                 ci.level = Omega)
  
  Opt.ests <-
    get_hpd_intervals_DossOpt(X = Samples,
                              h_level = Alpha,
                              ci_level = Omega)
  
  Batches.ests <- get_hpd_intervals_batches(X = Samples,
                                            h_level = Alpha,
                                            ci_level = Omega)
  
  concat <- do.call(
    rbind,
    list(
      Meeker.ests,
      SPIn.ests,
      Doss.ests,
      Doss.sbm.ests,
      Stan.ests,
      Opt.ests,
      Mixture.N.ests,
      Mixture.G.ests,
      Batches.ests
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

Mu <- 4.01
Sigma <- 0.0352
Type <- "none"
eff <- .02
Phi <- (1 - eff) / (eff + 1) ## if Phi = 0, independent samples
exper <- 3

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

Ncores <- 6
simu.time <- system.time(simus <- do.call(
  rbind,
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
                     }, mc.cores = Ncores)
))
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

ans$mu <- Mu
ans$sigma <- Sigma
ans$efficiency <- eff
ans$alpha <- Alpha
ans$ci_target_coverage <- TargetCoverage
ans$n_draws <- M
ans$asymmetry <- Type

ans

write.csv(ans,
          file = paste0("experiment_", exper, ".csv"),
          row.names = FALSE)

pp <- ggplot(simus.df,
             aes(x = point, fill = method)) +
  geom_density(alpha = 0.4) +
  facet_wrap(. ~ quantity, scales = "free") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(data = true.vals,
             aes(xintercept = true_value),
             linetype = "longdash") +
  theme_bw()

pp

ggsave(
  plot = pp,
  filename =
    paste0("../figures/HPD_experiment_",
           exper,
           ".pdf"),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)


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
  facet_grid(. ~ method) +
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
  facet_grid(. ~ method) +
  theme_bw(base_size = 20)

mean((simus$min < true.HPD[1]) * (simus$max > true.HPD[2]))
