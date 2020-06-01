\name{kott.quantile}
\alias{kott.quantile}
\title{Estimation of quantiles}
\description{
Calculates estimates, standard errors and confidence intervals for quantiles in subpopulations.
}
\usage{
kott.quantile(deskott, y, probs = c(0.25,0.50,0.75), by = NULL,
              vartype = c("se", "cv", "cvpct", "var"), 
              conf.int = FALSE, conf.lev = 0.95)
}
\arguments{
  \item{deskott}{Object of class \code{kott.design} containing the replicated survey data.}
  \item{y}{Formula defining the variable of interest.}
  \item{probs}{Vector of probability values to be used to calculate the quantiles estimates. The default value selects the quartiles estimates.} 
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL} (the default option) estimates refer to the whole population.}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose one or more of: standard error (the default), coefficient of variation, percent coefficient of variation, or variance.}
  \item{conf.int}{Boolean (\code{logical}) value to request confidence intervals for the estimates: the default is \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
}
\details{
This function calculates weighted estimates for the quantiles of a quantitative variable using suitable weights depending on the class of \code{deskott}: calibrated weights for class \code{kott.cal.design} and direct weights otherwise. Standard errors are calculated using the extended DAGJK method [Kott 99-01].

The mandatory argument \code{y} identifies the variable of interest, that is the variable for which quantiles estimates are to be calculated. The \code{deskott} variable referenced by \code{y} must be \code{numeric} and must not contain any missing value (\code{NA}).

The optional argument \code{probs} specifies the probability values (\code{0<=probs[i]<=1}) for which quantiles estimates must be calculated; the default option selects quartiles estimates. If \code{probs[i]} is equal to \code{0} (\code{1}) the corresponding "estimate" produced by \code{kott.quantile} coincides with the smallest (largest) observed value for the \code{y} variable.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{kottby} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. The \code{deskott} variables referenced by \code{by} (if any) must be \code{factor} and must not contain any missing value (\code{NA}).

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}. Given an input \code{kott.design} object with \code{nrg} random groups, \code{kott.quantile} builds the confidence intervals making use of a t distribution with \code{nrg-1} degrees of freedom.
}
\value{
The return value depends on the value of the input parameters. In the most general case, the function returns an object of class \code{list} (typically a list made up of data frames).
}
\section{Warning}{
It may happen that, in certain subpopulations, some of the \code{nrg} replicate weights turn out to be all zero: for these replicates it is not possible to provide quantiles estimates. In these cases, \code{kott.quantile} (i) returns \code{NaN} for the corresponding standard errors and (ii) prints a \code{warning} message.
}
\note{
Let \eqn{\hat{F}_y} be the estimate of the cumulative distribution of the \eqn{y} variabile. If an observed value \eqn{y^*} exists such that \eqn{\hat{F}_y(y^*)=probs[i]} than the i-th quantile estimate provided by \code{kott.quantile} equals \eqn{y^*}. If this is not the case, the \code{kott.quantile} function (i) finds the two observed values \eqn{y^-} and \eqn{y^+} (\eqn{y^- < y^+}) such that the corresponding values \eqn{\hat{F}_y(y^-)} and \eqn{\hat{F}_y(y^+)} are the closest to \eqn{probs[i]}, (ii) linearly interpolates \eqn{\hat{F}_y} between \eqn{\hat{F}_y(y^-)} and \eqn{\hat{F}_y(y^+)} and (iii) estimates the i-th quantile by inverting the linear approximation in the point \eqn{probs[i]}. 

The rigorous results of [kott 99-01] show that the DAGJK variance estimator for a given estimator \eqn{\hat{\theta}} is correct provided that PSUs are sampled with replacement and that \eqn{\hat{\theta}} is a smooth function of total estimators. As a result, it is not possible to guarantee that the DAGJK quantile variance estimator provided by \code{kott.quantile} is not biased.
}
\author{Diego Zardetto}
\references{
Kott, Phillip S. (1999) \emph{"The Extended Delete-A-Group Jackknife"}. Bulletin of the International Statistical Instititute. 52nd Session. Contributed Papers. Book 2, pp. 167-168.

Kott, Phillip S. (2001) \emph{"The Delete-A-Group Jackknife"}. Journal of Official Statistics, Vol.17, No.4, pp. 521-526.
}
\seealso{
\code{\link{kottby}} for estimating totals and means, \code{\link{kott.ratio}} for estimating ratios between totals, \code{\link{kott.regcoef}} for estimating regression coefficients and \code{\link{kottby.user}} for calculating estimates based on user-defined estimators.
}
\examples{
data(data.examples)

# Creation of a kott.design object:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)

# Estimate of the deciles of the income variable for
# the whole population:
kott.quantile(kdes,~income,probs=seq(0.1,0.9,0.1))

# Estimate of the median of income by age5c:
kott.quantile(kdes,~income,probs=0.5,by=~age5c,conf.int=TRUE)

# "Estimate" of the minimum and maximum of income by sex
# (notice the value of SE): 
kott.quantile(kdes,~income,probs=c(0,1),by=~sex)
}
\keyword{survey}