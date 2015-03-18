crisk.ncens.sim <- 
  function (foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0, 
            dist.ev, dist.cens, i, nsit) 
{
  nid    <- NA
  start  <- NA
  stop   <- NA
  obs    <- NA
  it     <- NA
  time   <- NA
  pro    <- vector()
  cause  <- NA
  a.cens <- NA
  b.cens <- NA
  obs[1] <- 1
  k.ev   <- 1
  sum    <- 0
  cshaz  <- list()
  az1    <- vector()
  
  if (is.null(z))
  {
    for (i in 1:nsit)
    {
      az1[i] <- 1
    }
  }else{
    for (i in 1:length(z))
    {
      if (!is.na(z[[i]][1]) && z[[i]][1] == "gamma") 
        az1[i] <- rgamma(1, as.numeric(z[[i]][2]), as.numeric(z[[i]][3]))
      if (!is.na(z[[i]][1]) && z[[i]][1] == "exp") 
        az1[i] <- rgamma(1, 1, as.numeric(z[[i]][2]))
      if (!is.na(z[[i]][1]) && z[[i]][1] == "weibull") 
        az1[i] <- rweibull(1, as.numeric(z[[i]][2]), as.numeric(z[[i]][3]))
      if (!is.na(z[[i]][1]) && z[[i]][1] == "unif") 
        az1[i] <- runif(1, as.numeric(z[[i]][2]), as.numeric(z[[i]][3]))
      if (!is.na(z[[i]][1]) && z[[i]][1] == "invgauss") 
        az1[i] <- rinvgauss(1, as.numeric(z[[i]][2]), as.numeric(z[[i]][3]))
    }
    if (length(z) == 1)
    {
      for (i in 2:nsit)
      {
        az1[i] <- az1[1]
      }
    }
  }
  
  if (dist.cens == "llogistic") {
    tc <- exp(rlogis(1, beta0.cens, anc.cens))
  }
  else {
    if (dist.cens == "weibull") {
      a.cens <- anc.cens
      b.cens <- (1/exp(-anc.cens * (beta0.cens)))^(1/anc.cens)
      tc <- rweibull(1, a.cens, b.cens)
    }
    else {
      if (dist.cens == "lnorm") {
        tc <- rlnorm(1, beta0.cens, anc.cens)
      }
    else {
      if (dist.cens== "unif") {
        tc <- runif(1, beta0.cens, anc.cens)
      }
    }
    }
  }
  suma <- 0
  if (!is.na(beta[1])) suma <- sum(sapply(beta, "[", seq(1,nsit,1)) * eff)
  # Cause-specific hazards
  for (k in 1:nsit)
  {
    if (dist.ev[k] == "llogistic") {
    cshaz[[k]] <- function(t) {return(az1[k]*exp(dlogis(t, beta0.ev[k] + suma, anc.ev[k])))}
    }
    else {
      if (dist.ev[k] == "weibull") {
        a.ev <- anc.ev
        b.ev <- (1/exp(-anc.ev * (beta0.ev+suma)))^(1/anc.ev)
        cshaz[[k]] <- function(t) {return(az1[k]*dweibull(t, a.ev[k], b.ev[k]))}
      }
      else {
        if (dist.ev[k] == "lnorm") {
          cshaz[[k]] <- function(t) {return(az1[k]*dlnorm(t, beta0.ev[k]+suma, anc.ev[k]))}
        }#if
      }#if
    }#if
  }#for
  
  A <- function(t,y){ #Cumulative all-cause hazard A
    suma <- 0
    for (k in 1:length(cshaz))
    {
      res <- suma + integrate(cshaz[[k]],lower=0,upper=t,subdivisions=1000)$value
    }
    res <- res + y
    return(res[1])
  }
  u     <- runif(1)
  iters <- 0
  while (A(0.0001,log(1-u))*A(foltime,log(1-u))>0 & iters < 1000)
  {
    u     <- runif(1)
    iters <- iters + 1
  }
  if (iters >= 1000) stop("Error: Values at endpoints not of opposite sign. \n")
  res <- uniroot(A, c(0.0001,foltime),tol=0.0001,y=log(1-u))
  tb  <- res$root
  
  sumprob <- 0
  for (k in 1:length(cshaz))
  {
    sumprob <- sumprob + cshaz[[k]](tb) 
  }
  for (k in 1:length(cshaz))
  {
      pro[k] <- cshaz[[k]](tb) / sumprob
  }
  cause1 <- rmultinom(1, 1, prob = pro)
  for (k in 1:length(cshaz))
  {
    if (cause1[k] == 1) cause <- k
  }
  az <- az1[k]
  nid <- i
  start <- 0
  it <- 0
  time <- tc
  if (tb < tc) {
      it <- 1
      time <- tb
    }
    stop <- time
    if (start < foltime && stop > foltime) {
      stop <- foltime
      time <- foltime
      it <- 0
    }
    if (start < 0 && stop > 0) {
      start <- 0
      time <- stop
    }
    
  sim.ind <- data.frame(nid = nid, cause = cause, time = time, status = it, 
                        start = start, stop = stop, z = az)
  for (k in 1:length(eff)) {
    sim.ind <- cbind(sim.ind, x = eff[k])
  }
  sim.ind <- subset(sim.ind, start < foltime & stop > 0)
  return(sim.ind)
}