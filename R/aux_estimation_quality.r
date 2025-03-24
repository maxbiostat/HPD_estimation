is_in <- function(x, l, u, finite = FALSE) {
  below <- x >= l
  above <- x <= u
  result <- as.logical(below * above)
  if (finite) {
    testL <- is.finite(l)
    testU <- is.finite(u)
    result <- as.logical(result * testL * testU)
  }
  return(result)
}
compute_MRAE <- function(hats, theta0) {
  mean(
    abs(hats - theta0)/abs(theta0)
  )
}
compute_MSE <- function(hats, theta0) {
  VV <- var(hats)
  BB <- mean(hats) - theta0
  MSE <- mean((hats - theta0) ^ 2)
  return(c(
    est_variance = VV,
    est_bias = BB,
    mse = MSE
  ))
}