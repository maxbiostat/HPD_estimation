library(tidyverse)

compute_coverage <- function(Mu,
                             Sigma,
                             eff,
                             M,
                             LvL = 0.95){
  
  estimates.df <- read_csv(
    file = paste0("../results/HPD_estimates_",
                  Mu,
                  "_s=", Sigma,
                  "_eff=", eff,
                  "_M=", M, "_bootEstimators.csv"))
  
  coverage <- aggregate(covers ~ quantity + type + nu,
                        subset(estimates.df), mean)
  coverage$coverage <- "full"
  
  
  coverage2 <- aggregate(covers ~ quantity + type + nu,
                         subset(estimates.df, min < true_value),
                         mean)
  coverage2$coverage <- "conditional"
  
  return(rbind(coverage, coverage2))
  
}


cc.extreme <- compute_coverage(Mu = 0, Sigma = 1,
                        eff = .05, M = 1e4,
                        LvL = .95)

cc.asym <- compute_coverage(Mu = 0.8314, Sigma = 0.30696,
                        eff = .05, M = 1e4,
                        LvL = .95)

cc.sym <- compute_coverage(Mu = 4.01, Sigma = 0.0345,
                            eff = .05, M = 1e4,
                            LvL = .95)


subset(cc.extreme, nu == "estimated" & coverage == "full")
subset(cc.asym, nu == "estimated" & coverage == "full")
subset(cc.sym, nu == "estimated" & coverage == "full")