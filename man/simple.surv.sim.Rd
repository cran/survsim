\encoding{utf8} 
\name{simple.surv.sim}
\alias{simple.surv.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a cohort with single-event survival times
}
\description{
Simulation of cohorts in a context of standard survival analysis including several covariates, individual heterogeneity and periods at risk prior and after the start of follow-up. 
}
\usage{
simple.surv.sim(n, foltime, dist.ev, anc.ev, beta0.ev, dist.cens="weibull", 
anc.cens, beta0.cens, z=NA, beta=NA, x=NA)
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
time to event distributions, with possible values \code{weibull} for the Weibull distribution, \code{lnorm} for the log-normal distribution and \code{llogistic} for the log-logistic distribution.
}
  \item{anc.ev}{
ancillary parameter for the time to event distribution.
}
  \item{beta0.ev}{
\eqn{\beta_0} parameter for the time to event distribution.
}
  \item{dist.cens}{
string indicating the time to censoring distribution, with possible values \code{weibull} for the Weibull distribution, \code{lnorm} for the log-normal distribution and \code{llogistic} for the log-logistic distribution. If no distribution is introduced, the time to censoring is assumed to follow a Weibull distribution.
}
  \item{anc.cens}{
real number containing the ancillary parameter for the time to censoring distribution.
}
  \item{beta0.cens}{
real number containing the \eqn{\beta_0} parameter for the time to censoring distribution.
}
  \item{z}{
vector with three elements that contains information relative to a random effect used in order to introduce individual heterogeneity. The first element indicates the distribution: \code{"unif"} states for a uniform distribution, \code{"gamma"} states for a gamma distribution, \code{"exp"} states for an exponential distribution, \code{"weibull"} states for a Weibull distribution and \code{"invgauss"} states for an inverse gaussian distribution. The second and third elements indicates the minimum and maximum in the case of a uniform distribution (both must be positive) and the parameters in the case of the rest of distributions. Notice that that just one parameter is needed in the case of exponential distribution. Its default value is \code{NA}, indicating that no individual heterogeneity is introduced.
}
  \item{beta}{
list of elements indicating the effect of the corresponding covariate. The number of vectors in \code{beta} must match the number of covariates. Its default value is NA, indicating that no covariates are included.
}
  \item{x}{
list of vectors indicating the distribution and parameters of any covariate that the user needs to introduce in the simulated cohort. The possible distributions are \code{"normal"} for normal distribution, \code{"unif"} for uniform distribution and \code{"bern"} for Bernoulli distribution. Its default value is NA, indicating that no covariates are included. The number of vectors in \code{x} must match the number of vectors in \code{beta}. Each vector in \code{x} must contain the name of the distribution and the parameter(s), which are: the probability of success in the case of a Bernoulli distribution; the mean and the variance in the case of a \code{normal} distribution; and the minimum and maximum in the case of a \code{uniform} distribution.
}

}

\value{
An object of class \code{mult.ev.data.sim}. It is a data frame containing the events suffered by the corresponding subjects. The columns of this data frame are detailed below
\item{nid }{an integer number that identifies the subject.}
\item{status }{logical value indicating if the corresponding event has been suffered or not.}
\item{start }{time at which the follow-up time begins for each event.}
\item{stop }{time at which the follow-up time ends for each event.}
\item{z}{Individual heterogeneity generated according to the specified distribution.}
\item{x}{value of each covariate randomly generated for each subject in the cohort.}
}
\references{
Kelly PJ, Lim LL. Survival analysis for recurrent event data: an application to childhood infectious diseases. Stat Med 2000 Jan 15;\bold{19}(1):13-33. 

}

\author{
David Moriña, Centre Tecnològic de Nutrició i Salut and Albert Navarro, Universitat Autònoma de Barcelona
}

\seealso{
\code{\link{survsim-package}}, \code{\link{accum}}, \code{\link{rec.ev.sim}}, \code{\link{mult.ev.sim}}
}
\examples{
### A cohort with 1000 subjects, with a maximum follow-up time of 3600 days and two 
### covariates, following a Bernoulli and uniform distribution respectively, and a 
### corresponding beta of -0.4 for the first covariate and a corresponding beta of 0
### for the second covariate. Notice that the time to censorship is assumed to 
### follow a Weibull distribution, as no other distribution is stated.

sim.data <- simple.surv.sim(n=1000, foltime=3600, dist.ev=c('llogistic'),
anc.ev=c(0.69978200185280),beta0.ev=c(5.84298525742252),,anc.cens=1.17783687569519,
beta0.cens=7.39773677281100,z=c("unif", 0.8,1.2), beta=list(c(-0.4),
c(0)), x=list(c("bern", 0.5), c("unif", 0.7, 1.3)))

summary(sim.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survsim }
\keyword{ individual }
\keyword{ single }
\keyword{ survival }
\keyword{ simulation }
