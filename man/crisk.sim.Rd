\encoding{utf8} 
\name{crisk.sim}
\alias{crisk.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a cohort in a competing risks context
}
\description{
Simulation of cohorts in a context of competing risks survival analysis including several covariates, individual heterogeneity and periods at risk prior and after the start of follow-up. 

Competing risks analysis considers time-to-first-event and the event type, possibly subject to right censoring (Beyersmann et al., 2009)

}
\usage{
crisk.sim(n, foltime, dist.ev, anc.ev, beta0.ev, dist.cens="weibull", 
anc.cens, beta0.cens, z=NULL, beta=NA, x=NA, nsit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
integer value indicating the desired size of the cohort to be simulated.
}
  \item{foltime}{
real number that indicates the maximum time of follow-up of the simulated cohort.
}
  \item{dist.ev}{
vector of arbitrary size indicating the time to event distributions, with possible values \code{weibull} for the Weibull distribution, \code{lnorm} for the log-normal distribution and \code{llogistic} for the log-logistic distribution.
}
  \item{anc.ev}{
vector of arbitrary size of real components containing the ancillary parameters for the time to event distributions.
}
  \item{beta0.ev}{
vector of arbitrary size of real components containing the \eqn{\beta_0} parameters for the time to event distributions.
}
  \item{dist.cens}{
string indicating the time to censoring distributions, with possible values \code{weibull} for the Weibull distribution, \code{lnorm} for the log-normal distribution, \code{llogistic} for the log-logistic distribution and \code{unif} for the uniform distribution. If no distribution is introduced, the time to censoring is assumed to follow a Weibull distribution.
}
  \item{anc.cens}{
real number containing the ancillary parameter for the time to censoring distribution or the maximum in case of uniform distributed time to censoring.
}
  \item{beta0.cens}{
real number containing the \eqn{\beta_0} parameter for the time to censoring distribution or the minimum in case of uniform distributed time to censoring.
}
  \item{z}{
list of vectors with three elements containing information relative to a random effect used in order to introduce individual heterogeneity. Each vector in the list refers to a possible competing risk, so the number of vectors must be equal to \code{nsit} or equal to 1 if the same random effect will be used for all the causes. The first element indicates the distribution: \code{"unif"} states for a uniform distribution, \code{"gamma"} states for a gamma distribution, \code{"exp"} states for an exponential distribution, \code{"weibull"} states for a Weibull distribution and \code{"invgauss"} states for an inverse gaussian distribution. The second and third elements indicates the minimum and maximum in the case of a uniform distribution (both must be positive) and the parameters in the case of the rest of distributions. Note that just one parameter is needed in the case of exponential distribution. Its default value is \code{NULL}, indicating that no individual heterogeneity is introduced. 
}
  \item{beta}{
list of vectors indicating the effect of the corresponding covariate. The number of vectors in \code{beta} must match the number of covariates, and the length of each vector must match the number of events considered. Its default value is NA, indicating that no covariates are included.
}
  \item{x}{
list of vectors indicating the distribution and parameters of any covariate that the user needs to introduce in the simulated cohort. The possible distributions are \code{"normal"} for normal distribution, \code{"unif"} for uniform distribution and \code{"bern"} for Bernoulli distribution. Its default value is NA, indicating that no covariates are included. The number of vectors in \code{x} must match the number of vectors in \code{beta}. Each vector in \code{x} must contain the name of the distribution and the parameter(s), which are: the probability of success in the case of a Bernoulli distribution; the mean and the variance in the case of a \code{normal} distribution; and the minimum and maximum in the case of a \code{uniform} distribution.
}

  \item{nsit}{
Number of different events that a subject can suffer. It must match the number of distributions specified in \code{dist.ev}.
}
}
\details{
In order to get the function to work properly, the length of the vectors containing the parameters of the time to event and the number of distributions indicated in the parameter \code{dist.ev} must be the same. 
}
\value{
An object of class \code{mult.ev.data.sim}. It is a data frame containing the events suffered by the corresponding subjects. The columns of this data frame are detailed below
\item{nid }{an integer number that identifies the subject.}
\item{cause }{cause of the event corresponding to the follow-up time of the individual.}
\item{time }{time until the corresponding event happens (or time to subject drop-out).}
\item{status }{logical value indicating if the corresponding event has been suffered or not.}
\item{start }{time at which the follow-up time begins for each event.}
\item{stop }{time at which the follow-up time ends for each event.}
\item{z }{Individual heterogeneity generated according to the specified distribution.}
\item{x }{value of each covariate randomly generated for each subject in the cohort.}
}
\references{
Beyersmann J, Latouche A, Buchholz A, Schumacher M. Simulating competing risks data in survival analysis. Stat Med 2009 Jan 5;\bold{28}(1):956-971. 

}

\author{
David Moriña, Centre for Research in Environmental Epidemiology (CREAL) and Albert Navarro, Universitat Autònoma de Barcelona
}

\seealso{
\code{\link{survsim-package}}, \code{\link{accum}}, \code{\link{rec.ev.sim}}, \code{\link{mult.ev.sim}}, \code{\link{simple.surv.sim}}
}
\examples{
### A cohort with 500 subjects, with a maximum follow-up time of 100 days and two 
### covariates, following Bernoulli distributions, and a corresponding beta of 
### 0.1698695 and 0.0007010932 for each event for the first covariate and a 
### corresponding beta of 0.3735146 and 0.5591244 for each event for the 
### second covariate. Notice that the time to censorship is assumed to follow a 
### log-normal distribution.

sim.data <- crisk.sim(n=500, foltime=100, dist.ev=c("lnorm","lnorm"),
anc.ev=c(1.479687, 0.5268302),beta0.ev=c(3.80342, 2.535374),dist.cens="lnorm",
anc.cens=1.242733,beta0.cens=5.421748,z=list(c("unif", 0.8,1.2), c("unif", 0.9, 1.5)), 
beta=list(c(0.1698695,0.0007010932),c(0.3735146,0.5591244)), 
x=list(c("bern", 0.381), c("bern", 0.564)), nsit=2)

summary(sim.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survsim }
\keyword{ individual }
\keyword{ competing risks }
\keyword{ survival }
\keyword{ simulation }
