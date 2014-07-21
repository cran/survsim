\encoding{utf8} 
\name{rec.ev.sim}
\alias{rec.ev.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a cohort with recurrent events
}
\description{
Simulation of cohorts in a context of recurrent event survival analysis including several covariates, individual heterogeneity and periods at risk before and after the initial time of follow-up.

Recurrent event data is a type of multiple event where the subject can experience repeated occurences of the same type (Kelly, 2000), for example repeated asthma attacks or sick leave episodes. In practice, the hazard of an recurrent event can vary depending on the number of previous occurrences, in terms of shape and intensity (Reis, 2011; Navarro, 2012). However, simulations based on a mixture of distributions with different baseline hazard rates are quite rare (Bender, 2005; Metcalfe, 2006).

In a recurrent data context, each subject can present different number of episodes. We talk of episodes (or occurrences) rather than events since each occurrence is a new episode of the same event. This package supposes that there exists one different and independent \eqn{Y_k} distribution depending on \eqn{k}, the number of episode at risk. The simulating process for each \eqn{Y_k} is the same than for the multiple events situation (see \code{\link{mult.ev.sim}}), but in this case, obviously, a subject cannot be at risk for the \eqn{k}-th episode if he/she hadn't had the \eqn{k-1}-th. 
}
\usage{
rec.ev.sim(n, foltime, dist.ev, anc.ev, beta0.ev, dist.cens=rep("weibull",
length(beta0.cens)), anc.cens, beta0.cens, z=NA, beta=NA, x=NA, lambda=NA, 
max.ep=Inf, priskb=0, max.old=0)
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
vector of arbitrary size indicating the time to event distributions, with possible values \code{weibull} for the Weibull distribution, \code{lnorm} for the log-normal distribution and \code{llogistic} for the log-logistic distribution. If a subject suffers more episodes than specified distributions, the last distribution specified here is used to generate times corresponding to posterior episodes.
}
  \item{anc.ev}{
vector of arbitrary size of real components containing the ancillary parameters for the time to event distributions.
}
  \item{beta0.ev}{
vector of arbitrary size of real components containing the \eqn{\beta_0} parameters for the time to event distributions.
}
  \item{dist.cens}{
vector of arbitrary size indicating the time to censoring distributions, with possible values \code{weibull} for the Weibull distribution, \code{lnorm} for the log-normal distribution and \code{llogistic} for the log-logistic distribution. If no distribution are introduced, the time to censoring is assumed to follow as many different Weibull distributions as different parameters are introduced in \code{beta0.cens}.
}
  \item{anc.cens}{
vector of arbitrary size of real components containing the ancillary parameters for the time to censoring distributions.
}
  \item{beta0.cens}{
vector of arbitrary size of real components containing the \eqn{\beta_0} parameters for the time to censoring distributions.
}
  \item{z}{
vector with three elements that contains information relative to a random effect used in order to introduce individual heterogeneity. The first element indicates the distribution: \code{"unif"} states for a uniform distribution, \code{"gamma"} states for a gamma distribution, \code{"exp"} states for an exponential distribution, \code{"weibull"} states for a Weibull distribution and \code{"invgauss"} states for an inverse gaussian distribution. The second and third elements indicates the minimum and maximum in the case of a uniform distribution (both must be positive) and the parameters in the case of the rest of distributions. Note that just one parameter is needed in the case of exponential distribution. Its default value is \code{NA}, indicating that no individual heterogeneity is introduced.
}
  \item{beta}{
list of vectors indicating the effect of the corresponding covariate. The number of vectors in \code{beta} must match the number of covariates, and the length of each vector must match the number of events considered. Its default value is NA, indicating that no covariates are included.
}
  \item{x}{
list of vectors indicating the distribution and parameters of any covariate that the user need to introduce in the simulated cohort. The possible distributions are \code{"normal"} for a normal distribution, \code{"unif"} for a uniform distribution and \code{"bern"} for a Bernoulli distribution. Its default value is NA, indicating that no covariates are included. The number of vectors in \code{x} must match the number of vectors in \code{beta}. Each vector in \code{x} must contain the name of the distribution and the parameter(s), which are: the probability of success in the case of a Bernoulli distribution; the mean and the variance in the case of a \code{normal} distribution; and the minimum and maximum in the case of a \code{uniform} distribution.
}
  \item{lambda}{
real number indicating the mean duration of each event or discontinous risk time, assumed to follow a zero-truncated Poisson distribution. Its default value is \code{NA}, corresponding to the case where the duration of each event or discontinous risk time is unnecessary information for the user.
}
  \item{max.ep}{
integer value that matches the maximum permitted number of episodes per subject. Its default value is \code{Inf}, i.e. the number of episodes per subject is no limited.
}
  \item{priskb}{
proportion of subjects at risk prior to the start of follow-up, defaults to 0.
}
  \item{max.old}{
maximum time at risk prior to the start of follow-up.
}
}
\details{
In order to get the function to work properly, the length of the vectors containing the parameters of the time to event and time to censure distributions and the number of distributions indicated in the parameter \code{dist} must be the same. Finally, \code{priskb} and \code{max.old} must be positive numbers, with \code{priskb} being between 0 and 1. Notice that large values of \code{max.old} can result in the routine taking a long time to simulate a cohort with the specified size.
}
\value{
An object of class \code{rec.ev.data.sim}. It is a data frame containing the episodes suffered by the corresponding subjects. The columns of the data frame are detailed below
\item{nid }{an integer number that identifies the subject.}
\item{real.episode }{number of the episode corresponding to the real history of the individual.}
\item{obs.episode }{number of the episode corresponding to the follow-up time of the individual.}
\item{time }{time until the corresponding event happens (or time to subject drop-out), regarding the beginning of the follow-up time.}
\item{status }{logical value indicating if the episode corresponds to an event or a drop-out.}
\item{start }{time at which an episode starts, taking the beginning of follow-up as the origin of the time scale.}
\item{stop }{time at which an episode ends, taking the beginning of follow-up as the origin of the time scale.}
\item{time2 }{time until the corresponding event happens (or time to subject drop-out), in calendar time.}
\item{start2 }{time at which an episode starts, where the time scale is calendar time.}
\item{stop2 }{time at which an episode ends, where the time scale is calendar time.}
\item{old }{real value indicating the time that the individual was at risk before the beginning of follow-up.}
\item{risk.bef }{factor that indicates if an individual was at risk before the beginning of follow-up or not.}
\item{long }{time not at risk immediately after an episode.}
\item{z}{Individual heterogeneity generated according to the specified distribution.}
\item{x}{value of each covariate randomly generated for each subject in the cohort.}
}
\references{
Kelly PJ, Lim LL. Survival analysis for recurrent event data: an application to childhood infectious diseases. Stat Med 2000 Jan 15;\bold{19}(1):13-33. 

Bender R, Augustin T, Blettner M. Generating survival times to simulate Cox proportional hazards models. Stat Med 2005 Jun 15;\bold{24}(11):1713-1723.

Metcalfe C, Thompson SG. The importance of varying the event generation process in simulation studies of statistical methods for recurrent events. Stat Med 2006 Jan 15;\bold{25}(1):165-179.

Reis RJ, Utzet M, La Rocca PF, Nedel FB, Martin M, Navarro A. Previous sick leaves as predictor of subsequent ones. Int Arch Occup Environ Health 2011 Jun;\bold{84}(5):491-499.

Navarro A, Moriña D, Reis R, Nedel FB, Martin M, Alvarado S. Hazard functions to describe patterns of new and recurrent sick leave episodes for different diagnoses. Scand J Work Environ Health 2012 Jan 27.

Moriña D, Navarro A. The R package survsim for the simulation of simple and complex survival data. Journal of Statistical Software 2014 Jul; \bold{59}(2):1-21.

}

\author{
David Moriña, Centre for Research in Environmental Epidemiology (CREAL) and Albert Navarro, Universitat Autònoma de Barcelona
}

\seealso{
\code{\link{survsim-package}}, \code{\link{accum}}, \code{\link{mult.ev.sim}}, \code{\link{simple.surv.sim}}, \code{\link{crisk.sim}}
}
\examples{
### A cohort with 500 subjects, with a maximum follow-up time of 1825 days and 
### just a covariate, following a Bernoulli distribution, and a corresponding 
### beta of -0.4, -0.5, -0.6 and -0.7 for each episode.

sim.data <- rec.ev.sim(n=500, foltime=1825, dist.ev=c('lnorm','llogistic', 
'weibull','weibull'),anc.ev=c(1.498, 0.924, 0.923, 1.051),beta0.ev=c(7.195, 
6.583, 6.678, 6.430),,anc.cens=c(1.272, 1.218, 1.341, 1.484),
beta0.cens=c(7.315, 6.975, 6.712, 6.399), z=c("unif", 0.8,1.2), 
beta=list(c(-0.4,-0.5,-0.6,-0.7)), x=list(c("bern", 0.5)),
lambda=c(2.18,2.33,2.40,3.46), priskb=0.5, max.old=730)

summary(sim.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survsim }
\keyword{ individual }
\keyword{ recurrent }
\keyword{ survival }
\keyword{ simulation }
