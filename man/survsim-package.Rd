\encoding{utf8} 
\name{survsim-package}
\alias{survsim-package}
\alias{survsim}
\docType{package}
\title{
Simulation of simple and complex survival data 
}
\description{
Simulation of cohorts in a context of simple and complex survival analysis, multiple events and recurrent events including several covariates, individual heterogeneity and periods at risk before and after the initial time of follow-up.

\tabular{cccc}{
         Distribution \tab Survival function   \tab Density function \tab Parametrization\cr
         Weibull \tab \eqn{exp(- \lambda t^p)} \tab \eqn{\lambda pt^{p-1}exp(- \lambda t^p)} \tab \eqn{\lambda = exp(-p \beta_0)}\cr
         Log-normal \tab \eqn{1- \Theta((log(t)- \mu)/ \sigma)} \tab \eqn{(1/(t \sigma \sqrt{2 \pi})) exp((-1/(2 \sigma^2))(log(t) - \mu)^2)} \tab \eqn{\mu = \beta_0}\cr
         Log-logistic \tab \eqn{1/(1+(\lambda t)^{1/ \gamma}}) \tab \eqn{\lambda^{1/ \gamma}t^{1/ (\gamma - 1)}/ (\gamma (1 + (\lambda t)^{1/ \gamma})^2)} \tab \eqn{\lambda = exp(- \beta_0)}\cr
       }
\tabular{cc}{
         Distribution \tab Time\cr
         Weibull \tab \eqn{t = (- ln u/ \lambda)^{1/p}}\cr
         Log-normal \tab \eqn{t = exp(\beta_0 + \gamma (log(u)-log(1-u)))}\cr
         Log-logistic \tab \eqn{t = exp(\beta_0 + \sigma \Theta^{-1}(u))}\cr
       }
       
Where \eqn{\Theta} is the standard normal cumulative distribution.

In order to simulate censored survival data, two survival distributions are required, one for the uncensored survival times that would be observed if the follow-up had been sufficiently long to reach the event and another representing the censoring mechanism. The uncensored survival distribution, \eqn{T'_i}, for \eqn{i=1,\ldots,n} subjects, could be generated to depend on a set of covariates with a specified relationship with survival, which represents the true prognostic importance of each covariate (Burton, 2006). The package allows to simulate times by means of using Weibull (and exponential as a particular case), log-normal and log-logistic distributions, as such is showed in previous table.
To induce individual heterogeneity or within-subject correlation we generate \eqn{Z_i}, a random effect covariate that follows a particular distribution (Uniform or Normal). 

\deqn{t_i = t_i'\cdot z_i}

When \eqn{z_i = 1}, for all subjects, we are in the case of individual homogeneity and the survival times are completely specified by the covariates.
Random non-informative right censoring, \eqn{C_i}, can be generated in a similar manner to the uncensored survival times, \eqn{T'_i}, by assuming a  particular distribution for the censoring times (previous table), but without including any covariates nor individual heterogenity. 
The observation times, \eqn{Y_i'}, incorporating both events and censored observations are calculated for each case by combining the uncensored survival times, \eqn{T_i}, and the censoring times, \eqn{C_i}. If the uncensored survival time for an observation is less than or equal to the censored time, then the event is considered to be observed and the observation time equals the uncensored survival time, otherwise the event is considered censored and the observation time equals the censored time. In other words, once simulated \eqn{t_i} and \eqn{c_i}, we can define \eqn{Y_i'= min(t_i,c_i)} as the obervation time with \eqn{\delta_i} an indicator of non-censoring, i.e. \eqn{\delta_i = I(t_i \le c_i )}.
While all \eqn{y_i'}  start at 0, the package allows create dynamic cohorts. We can generate entry times higher than 0 adding a \eqn{t_0} value corresponding with an uniform distribution in \eqn{[0,t_{max follow-up}]}. We can also simulate subjects at risk before of the initial time of follow-up \eqn{(y_i'= 0)}, by including an uniform distribution for \eqn{t_0} between \eqn{[-t_{max old},0)} for a fixed percentage of subjects. Then:

\deqn{y_i=y_i' + t_0}

where \eqn{t_0} follows a uniform distribution in \eqn{[0,t_{max follow-up}]} if entry time is 0 or more and \eqn{t_0} is uniform distributed in \eqn{[-t_{max old}, 0)} if entry time is less than 0.
Therefore, \eqn{t_0} represents the initial point of the episode, \eqn{y_i} the endpoint and \eqn{y_i'} is the lenght. Note that \eqn{y_i'+t_0} can be higher than \eqn{t_{max follow-up}}, and in this case \eqn{y_i} will be set at \eqn{t_{max follow-up}} and \eqn{\delta_i=0}. The observations corresponding to the subjects at risk before of the initial time of follow-up have \eqn{t_0} negative, then the initial point of the episode will be set at 0. \eqn{y_i} may also be negative, in this case this episode will not be included in the simulated data, as long as this episode won't be observed in practice. 

}
\details{
\tabular{ll}{
Package: \tab survsim\cr
Type: \tab Package\cr
Version: \tab 1.1.5\cr
Date: \tab 2018-05-17\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
The package provide a tool for simulation of cohorts in a simple single-event context through the function \code{\link{simple.surv.sim}}, in a recurrent event context with the function \code{\link{rec.ev.sim}}, in a multiple event context with the function \code{\link{mult.ev.sim}} and in a competing risks context with the function \code{\link{crisk.sim}}, and it also allows the user to generate aggregated data from the simulated cohort, by means of the function \code{\link{accum}}.
}
\author{
David Moriña, Centre for Research in Environmental Epidemiology (CREAL) and Albert Navarro, Universitat Autònoma de Barcelona

Maintainer: David Moriña <david.morina@uab.cat>
}

\references{

Kelly PJ, Lim LL. Survival analysis for recurrent event data: an application to childhood infectious diseases. Stat Med 2000 Jan 15;\bold{19}(1):13-33.

Bender R, Augustin T, Blettner M. Generating survival times to simulate Cox proportional hazards models. Stat Med 2005 Jun 15;\bold{24}(11):1713-1723.

Metcalfe C, Thompson SG. The importance of varying the event generation process in simulation studies of statistical methods for recurrent events. Stat Med 2006 Jan 15;\bold{25}(1):165-179.

Burton A, Altman DG, Royston P, Holder RL. The design of simulation studies in medical statistics. Stat Med 2006 Dec 30;\bold{25}(24):4279-4292.

Beyersmann J, Latouche A, Buchholz A, Schumacher M. Simulating competing risks data in survival analysis. Stat Med 2009 Jan 5;\bold{28}(1):956-971. 

Reis RJ, Utzet M, La Rocca PF, Nedel FB, Martin M, Navarro A. Previous sick leaves as predictor of subsequent ones. Int Arch Occup Environ Health 2011 Jun;\bold{84}(5):491-499.

Navarro A, Moriña D, Reis R, Nedel FB, Martin M, Alvarado S. Hazard functions to describe patterns of new and recurrent sick leave episodes for different diagnoses. Scand J Work Environ Health 2012 Jan 27.

Moriña D, Navarro A. The R package survsim for the simulation of simple and complex survival data. Journal of Statistical Software 2014 Jul; \bold{59}(2):1-20.

}
\keyword{ survsim }

