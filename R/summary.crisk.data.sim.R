summary.crisk.data.sim <-
function(object, ...)
  {
    if(!inherits(object, "crisk.data.sim")) stop("Wrong data type")
    sub.risk    <- vector()
    num.events  <- vector()
    foltime     <- vector()
    med.foltime <- vector()
    mean.ep.sub <- vector()
    dens.incid  <- vector()
    cause       <- vector()
    for (i in 1:attr(object,"nsit"))
    {
       cause[i]       <- i
       sub.risk[i]    <- attr(object,"n")
       num.events[i]  <- as.integer(sum(!is.na(object$status[object$cause==i])))
       foltime[i]     <- sum(object$time[object$cause==i & !is.na(object$cause)])
       med.foltime[i] <- median(object$time[object$cause==i & !is.na(object$cause)])
       mean.ep.sub[i] <- sum(object$status[object$cause==i & !is.na(object$cause)])/dim(object)[1]
       dens.incid[i]  <- num.events[i]/foltime[i]
    }
    ans <- data.frame(cause, sub.risk, num.events, foltime, med.foltime,
                      mean.ep.sub, dens.incid)
    class(ans)  <- "summary.crisk.data.sim"
    return(ans)
  }
