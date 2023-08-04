library(tidyverse)
library(devtools)
source_url("https://raw.githubusercontent.com/maxbiostat/ESS_HPDs/main/R/aux.r")
source_url("https://raw.githubusercontent.com/maxbiostat/ESS_HPDs/main/R/aux_fake_MCMC.r")
source("aux_estimation_quality.r")
##############################
############# Functions

est_icdf <- function(p, samples,
                     alpha = .95, method = "bm") {
  
  Z <- qnorm(p = (1 + alpha) / 2)
  
  ests <- mcmcse::mcse.q(x = samples,
                         q = p, method = method)
  
  CI <- ests$est[1] + Z * c(-1, 1) * ests$se[1]
  
  out <- tibble::tibble(
    p = p,
    xi_hat = ests$est,
    lwr = CI[1],
    upr = CI[2],
    method = method
  )
  
  return(out)
}

estimated_inverse_cdf_Doss <- function(x, ps,
                                       alpha = 0.95, method = "bm") {
  
  res <- do.call(rbind,
                 lapply(ps, function(p) {
                   est_icdf(p = p,
                            samples = x,
                            alpha = alpha,
                            method = method)
                 }))
  
  return(res)
  
}

run_once <- function(i, Mu, Sigma, Phi, Alpha, M, ps, logn) {
  
  Samples <- generate_fake_MCMC(
    N = M,
    mu = Mu,
    sigma = Sigma,
    phi = Phi,
    LN = logn
  )
  
  all.ests <- Reduce(rbind,
                 list(
                   estimated_inverse_cdf_Doss(x = Samples,
                                              ps = ps,
                                              alpha = Alpha, method = "bm"),
                   estimated_inverse_cdf_Doss(x = Samples,
                                              ps = ps,
                                              alpha = Alpha, method = "obm"),
                   estimated_inverse_cdf_Doss(x = Samples,
                                              ps = ps,
                                              alpha = Alpha, method = "sub")
                 ))
  out <- tibble::tibble(
    min = min(Samples),
    max = max(Samples),
    all.ests,
    replicate = i
  )
  return(out)
}
##########################################
### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
## should the target be log-normal (TRUE) or normal (FALSE)?
logNormal <- TRUE

Mu <- 0
Sigma <- 1
eff <- .02
Phi <- (1 - eff) / (1 + eff) ## if Phi = 0, independent samples

Alpha <- 0.95
TargetCoverage <- 0.95
M <- 1E6
Nrep <- 500

HPD.probs <- lognormal_hpd_percentiles(alpha = .95,
                                       lmean = Mu, lsd = Sigma)
pp <- c(HPD.probs[1],
        .025, .05, .1, .25, .5, .75, .90, .95, .975)

Ncores <- 10
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
                         ps = pp,
                         logn = logNormal
                       )
                     }, mc.cores = Ncores)
))
simu.time

true.vals <- tibble::tibble(p = pp,
                            true_value = qlnorm(p = pp, 
                                                meanlog = Mu,
                                                sdlog = Sigma))

simus.df <- merge(simus, true.vals, by = "p")
simus.df <-
  simus.df %>% mutate(covers = is_in(true_value, lwr, upr))
simus.df <-
  simus.df %>% mutate(lwr = ifelse(is.infinite(lwr), 0, lwr))

coverage <-
  aggregate(covers ~ p + method, simus.df, mean) ## coverage

ci.width <-
  aggregate((upr - lwr)/true_value ~ p + method, simus.df, mean) ## average CI width

bias <-
  aggregate((xi_hat - true_value) ~ p + method, simus.df, mean)

variance <- aggregate(xi_hat ~ p + method, simus.df, var)
mse <-
  aggregate((xi_hat - true_value) ^ 2 ~ p + method, simus.df, mean)

ans <-
  Reduce(function(...)
    merge(..., by = c("p", "method")),
    list(coverage, ci.width,
         bias,
         variance, mse))
names(ans) <- c("p",
                "method",
                "est_coverage",
                "ci_width",
                "bias",
                "variance",
                "MSE")

ans

write.csv(ans,
          paste0("saved_data/bigM_quantile_estimation_v4_Doss_BM_m=", Mu,
                 "_s=_", Sigma, "_eff=", eff, ".csv" ),
          row.names = FALSE)

library(scales)

covers <- ggplot(ans,
                 aes(x = p, y = est_coverage)) +
  geom_point() +
  # scale_x_log10(expression(p),
  #                    breaks = trans_breaks("log10", function(x) 10^x),
  #                    labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(expression(p)) + 
  scale_y_continuous("Attained coverage") + 
  geom_hline(yintercept = Alpha, linetype = "longdash") +
  facet_grid(method~.) + 
  ggtitle("Coverage") +
  theme_bw(base_size = 20)

covers

widths <- ggplot(ans,
                 aes(x = p, y = ci_width)) +
  geom_point() +
  # scale_x_log10(expression(p),
  #                    breaks = trans_breaks("log10", function(x) 10^x),
  #                    labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(expression(p)) + 
  scale_y_continuous("Confidence interval width") + 
  facet_grid(method~.) + 
  ggtitle("CI witdth") +
  theme_bw(base_size = 20)

widths

ggsave(
  plot = covers,
  filename = paste0("figures/bigM_quantile_coverage_v4_Doss_BM_m=", Mu,
                    "_s=_", Sigma, "_eff=", eff, ".pdf" ),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)

ggsave(
  plot = covers,
  filename = paste0("figures/bigM_quantile_coverage_v4_Doss_BM_m=", Mu,
                    "_s=_", Sigma, "_eff=", eff, ".pdf" ),
  scale = 1,
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)
