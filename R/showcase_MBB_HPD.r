library(boot)

Mu <-  0#0.8314 # 4.01 
Sigma <- 1# 0.30696  # 0.0345
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

j <- 234
samples <- Simus[[j]][1:M]
M <- length(samples)
B <- 500
nu <- 1/2
bsize <- round(M^nu)

HPD.fun.2 <- function(tsb) {
  c(HDInterval::hdi(tsb, credMass = LvL))
}

hpd.hat <- HPD.fun.2(samples)
rg <- round(range(samples), 3)

bb <- boot::tsboot(tseries = samples,
                   statistic = HPD.fun.2,
                   R = B,
                   l = bsize,
                   sim = "fixed",
                   parallel = "multicore", ncpus = 8)
bb
true.HPD
rg 
( a.CI <- boot::boot.ci(bb, index = 1,
                        type = c("basic", "norm", "perc"),
                        conf = TargetCoverage) )
( b.CI <- boot::boot.ci(bb, index = 2,
                        type = c("basic", "norm", "perc"),
                        conf = TargetCoverage) )


plot(bb$t, 
     xlab = expression(hat(a)),
     ylab = expression(hat(b)),
     main = paste0("Range = (", rg[1], ", ", rg[2], ")"),
     cex = .7)
points(x = hpd.hat[1], y = hpd.hat[2], col = 2, pch = 16, cex = 2)
points(x = true.HPD[1], y = true.HPD[2], col = "purple", pch = 18, cex = 2)
legend(x = "bottomright",
       pch = c(1, 16, 18),
       col = c("black", "red", "purple"),
       legend  = c("Bootstrap resamples",
                   "Observed",
                   "True"),
       bty = 'n')

apply(bb$t, 2, hist)