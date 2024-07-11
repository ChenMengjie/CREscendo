to_probs <- function (x) {
  n <- sum(x)
  return(x/n)
}

chisq_stat  <- function (obs, exp) {
  pr <- to_probs(exp)
  nexp <- pr * sum(obs)
  return(sum((obs - nexp)^2/nexp))
}

### output the chisq statistics for each group of sites
contigency <- function(obs, exp) {
  pr <- to_probs(exp)
  nexp <- pr * sum(obs)
  return((obs - nexp)/sqrt(nexp))
}

