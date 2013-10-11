rec.ev.sim <-
function(n, foltime, dist.ev, anc.ev, beta0.ev, dist.cens=rep("weibull",length(beta0.cens)), anc.cens, beta0.cens, z=NA, beta=NA, x=NA,
                        lambda=NA, max.ep=Inf, priskb=0, max.old=0)
{
  # Arguments check
  if (length(anc.ev) != length(beta0.ev)) stop("Wrong number of parameters")
  if (length(anc.cens) != length(beta0.cens)) stop("Wrong number of parameters")
  if (length(anc.ev) != length(dist.ev)) stop("Wrong number of parameters")
  if (priskb > 1 || priskb < 0) stop("Wrong proportion of left-censured individuals")
  if (max.old < 0) stop("Wrong maximum time before follow-up")
  if (!is.na(z) && length(z) != 3 && z[1] != "exp") stop("Wrong numbers of elements in z")
  if (!is.na(z) && length(z) != 2 && z[1] == "exp") stop("Wrong numbers of elements in z")
  if (!is.na(z) && z[1] == "unif" && as.numeric(z[3]) <= as.numeric(z[2])) stop("Wrong specification of z")
  if (!is.na(z) && z[1] == "unif" && as.numeric(z[2]) < 0) stop("Wrong specification of z")
  if (!is.na(x) && is.na(beta)) stop("Wrong specification of covariables!")
  if (is.na(x) && !is.na(beta)) stop("Wrong specification of covariables")
  
  sim.data <- list()
  eff      <- vector()
  un.ncens <- runif(n, 0, foltime)
  un.cens  <- runif(n, -max.old, 0)
  
  ncens    <- as.integer(priskb*n)
  nncens   <- as.integer(n - ncens)
  max.time <- max(foltime, max.old)
  
  eff[1]   <- 0
   
  if (nncens != 0)
  {
    for (i in 1:nncens)
    {
        if (!is.na(x[1]))
        {
          for (k in 1:length(x))
          {
            if (x[[k]][1] == "unif")   eff[k] <- runif(1,as.numeric(x[[k]][2]),as.numeric(x[[k]][3]))
            if (x[[k]][1] == "normal") eff[k] <- rnorm(1,as.numeric(x[[k]][2]),as.numeric(x[[k]][3]))
            if (x[[k]][1] == "bern")   eff[k] <- rbinom(1,1,as.numeric(x[[k]][2]))
          } #for
      } #if
      sim.data[[i]] <- rec.ev.ncens.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z, beta, eff, 
                                         lambda, dist.ev, dist.cens, max.ep, un.ncens[i], i, max.time)
    } #for
  } # if
  if (ncens != 0)
  {
    i <- nncens + 1
    repeat
    {
      if (!is.na(x[1]))
      {
        for (k in 1:length(x))
        {
          if (x[[k]][1] == "unif")   eff[k] <- runif(1,as.numeric(x[[k]][2]),as.numeric(x[[k]][3]))
          if (x[[k]][1] == "normal") eff[k] <- rnorm(1,as.numeric(x[[k]][2]),as.numeric(x[[k]][3]))
          if (x[[k]][1] == "bern")   eff[k] <- rbinom(1,1,as.numeric(x[[k]][2]))
        } #for
      } #if
      sim.data[[i]] <- rec.ev.cens.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z, beta, eff,
                                        lambda, dist.ev, dist.cens, max.ep, un.cens[i],i, max.time)
      if (dim(sim.data[[i]])[1] != 0)
      {
        i <- i + 1
      }
      
      if (i == n + 1) break
      
    } #repeat
  } #if
  sim.data <- do.call(rbind, sim.data)
  class(sim.data) <- c("rec.ev.data.sim", "data.frame")
  attr(sim.data, "n") <- n
  attr(sim.data, "foltime") <- foltime
  attr(sim.data, "ndist") <- length(dist.ev)
  return(sim.data)
}
