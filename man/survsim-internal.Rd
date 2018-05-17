\encoding{utf8} 
\name{survsim-internal}
\alias{survsim-internal}
\alias{rtpois}
\alias{simple.ev.sim}
\alias{crisk.ncens.sim}
\alias{mult.ev.ncens.sim}
\alias{rec.ev.cens.sim}
\alias{rec.ev.ncens.sim}
\alias{summary.simple.surv.sim}
\alias{print.summary.simple.surv.sim}
\alias{summary.crisk.data.sim}
\alias{print.summary.crisk.data.sim}
\alias{summary.mult.ev.data.sim}
\alias{print.summary.mult.ev.data.sim}
\alias{summary.rec.ev.data.sim}
\alias{print.summary.rec.ev.data.sim}

\title{Internal survsim functions}
\description{Internal survsim functions}
\usage{
rtpois(N, lambda)
simple.ev.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0,
                  dist.ev, dist.cens, i)
crisk.ncens.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0,
                  dist.ev, dist.cens, i, nsit)
mult.ev.ncens.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0,
                  dist.ev, dist.cens, i, nsit)
rec.ev.cens.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0,
                lambda=NA, dist.ev, dist.cens, max.ep=Inf, un.cens, i, max.time)
rec.ev.ncens.sim(foltime, anc.ev, beta0.ev, anc.cens, beta0.cens, z=NULL, beta=0, eff=0, 
                 lambda=NA, dist.ev, dist.cens, max.ep=Inf, un.ncens, i, max.time)
\method{summary}{simple.surv.sim}(object, ...)
\method{summary}{mult.ev.data.sim}(object, ...)
\method{summary}{rec.ev.data.sim}(object, ...)
\method{summary}{crisk.data.sim}(object, ...)
\method{print.summary}{simple.surv.sim}(x, ...)
\method{print.summary}{mult.ev.data.sim}(x, ...)
\method{print.summary}{rec.ev.data.sim}(x, ...)
\method{print.summary}{crisk.data.sim}(x, ...)
}

\details{These functions are not to be called by the user}
\seealso{
\code{\link{survsim-package}}, \code{\link{accum}}
}
\keyword{internal}
