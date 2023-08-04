### Asymmetric : (2.58, 0.254)
### Symmetric : (4.01, 0.0352)
Mu <- 0
Sigma <- 1
Alpha <- .95
Level <- .95

true.HPD <- lognormal_hpd(alpha = Alpha, lmean = Mu, lsd = Sigma)
true.Hprobs <- plnorm(true.HPD, m = Mu, s = Sigma)
true.ps <- c( c(1-Alpha, 1 + Alpha)/2, true.Hprobs, 1/2)  
ESSs <- c(50, 100, 200, 625, 1000)