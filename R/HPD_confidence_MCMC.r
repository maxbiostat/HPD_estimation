library(tidyverse)
library(fakeMCMC)
library(boot)

source("aux_estimation_quality.r")
source("aux_HPD_estimators.r") ## basically old estimators
##############################
############# Functions
run_once <- function(i, Mu, Sigma, LvL, Omega, M) {
  
  Samples <- Simus[[i]]
  
  concat <- do.call(rbind,
                    list(
                         hpd_intervals_boot(samples = Samples,
                                            alpha = LvL, gamma = TargetCoverage,
                                            B = 300, cores = 8),
                         hpd_intervals_boot(samples = Samples,
                                            alpha = LvL, gamma = TargetCoverage,
                                            B = 300, nu = 1/2,
                                            cores = 8),
                         hpd_intervals_boot(samples = Samples,
                                            alpha = LvL, gamma = TargetCoverage,
                                            B = 300, nu = 1/5,
                                            cores = 8)
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

Mu <-  0.8314 # 4.01 
Sigma <- 0.30696 # 0.0345
eff <- .05
LvL <- 0.95 ## HPD level
TargetCoverage <- 0.95 ## confidence interval
M <- 1E4
Nrep <- 500

true.HPD <- as.numeric(
  HDInterval::hdi(function(p) qlnorm(p, Mu, Sigma), credMass = LvL)
)
true.qs <- plnorm(q = true.HPD, Mu, Sigma)

load(
  paste0("../saved_data/LAR_simus_m=",
         Mu,
         "_s=_", Sigma,
         "_eff=", eff,
         "_M=", M,
         ".RData" )
)

Ncores <- 10
est.time <- system.time(estimates <- do.call(
  rbind,
  parallel::mclapply(seq_len(Nrep),
                     function(i) {
                       run_once(
                         i,
                         M = M,
                         Mu = Mu,
                         Sigma = Sigma,
                         LvL = LvL,
                         Omega = TargetCoverage
                       )
                     }, mc.cores = Ncores)
))
est.time

true.vals <- tibble::tibble(true_value = true.HPD,
                            quantity = c("HPD_L", "HPD_U"))

estimates.df <- merge(estimates, true.vals, by = "quantity")
estimates.df <-
  estimates.df %>% mutate(covers = is_in(true_value, lwr, upr))
estimates.df <-
  estimates.df %>% mutate(lwr = ifelse(is.infinite(lwr), 0, lwr))

coverage <- aggregate(covers ~ quantity + type + nu,
                      subset(estimates.df), mean)

subset(coverage, quantity == "HPD_L")
subset(coverage, quantity == "HPD_U")

coverage2 <- aggregate(covers ~ quantity + type + nu,
                       subset(estimates.df, min < true_value), mean)

subset(coverage2, quantity == "HPD_L")
subset(coverage2, quantity == "HPD_U")

head(estimates.df)
head(subset(estimates.df,
            type == "percent"))
head(subset(estimates.df,
            type == "percent" & min < true_value))

write.csv(estimates.df,
          file = paste0("../results/HPD_estimates_",
                        Mu,
                        "_s=", Sigma,
                        "_eff=", eff,
                        "_M=", M, "_bootEstimators.csv"),
          row.names = FALSE)