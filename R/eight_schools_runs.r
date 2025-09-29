library(cmdstanr)

golden <- cmdstanr_example(
  example = c("schools_ncp"),
  method = c("sample"),
  iter_sampling = 1e6,
  quiet = TRUE
)

bf.tau <- golden$draws("tau")

hist(bf.tau, probability = TRUE,
     main = "",
     xlab = expression(tau))

Alphas <- c(.5, .8, .9, .95)

true.HPD <- do.call(rbind,
                    lapply(Alphas, function(a){
                      hpd <- HDInterval::hdi(bf.tau, credMass = a)
                      out <- data.frame(HPD_L = hpd[1],
                                        HPD_U = hpd[2],
                                        alpha = a)
                      return(out)
                    }))
rownames(true.HPD) <- NULL

Fhat <- ecdf(bf.tau)

ps <- data.frame(t(apply(true.HPD, 1, function(x) Fhat(x[1:2]))))
colnames(ps) <- c("p_L", "p_U")

write.csv(cbind(true.HPD, ps),
          file = "../results/true_HPD_EightSchools.csv",
          row.names = FALSE)


generate_tau_draws <- function(M){
  out <- cmdstanr_example(
    example = c("schools_ncp"),
    method = c("sample"),
    chains = 1,
    iter_sampling = M,
    quiet = TRUE
  )
  return(out$draws("tau"))
}

Nrep <- 500
Simus <- vector(Nrep, mode = "list")
M <- 1e4

for(i in 1:Nrep){
  Simus[[i]] <- generate_tau_draws(M = M)
}


save(Simus,
     file = paste0("../saved_data/EightSchools_simus",
                   "_M=", M,
                   ".RData" ))
