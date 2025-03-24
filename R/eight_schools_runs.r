library(cmdstanr)

golden <- cmdstanr_example(
  example = c("schools_ncp"),
  method = c("sample"),
  iter_sampling = 1e6,
  quiet = TRUE
)

bf.tau <- golden$draws("tau")
hist(bf.tau)
( brute.force.HPD <- HDInterval::hdi(bf.tau) )

Fhat <- ecdf(tau)

Fhat(brute.force.HPD)

generate_tau_draws <- function(M){
  out <- cmdstanr_example(
    example = c("schools_ncp"),
    method = c("sample"),
    iter_sampling = M/4,
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
     file = paste0("saved_data/EightSchools_simus",
                   "_M=", M,
                   ".RData" ))
