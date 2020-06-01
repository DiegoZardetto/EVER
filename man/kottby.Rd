\name{kottby}
\alias{kottby}
\title{Estimation of totals and means}
\description{
Calculates estimates, standard errors and confidence intervals for totals and means in subpopulations.
}
\usage{
kottby(deskott, y, by = NULL, estimator = c("total", "mean"),
       vartype = c("se", "cv", "cvpct", "var"), 
       conf.int = FALSE, conf.lev = 0.95)
}
\arguments{
  \item{deskott}{Object of class \code{kott.design} containing the replicated survey data.}
  \item{y}{Formula defining the variables of interest.}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL} (the default option) estimates refer to the whole population.}
  \item{estimator}{\code{character} specifying the desired estimator: it may be \code{"total"} (the default) or \code{"mean"}.}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose one or more of: standard error (the default), coefficient of variation, percent coefficient of variation, or variance.}
  \item{conf.int}{Boolean (\code{logical}) value to request confidence intervals for the estimates: the default is \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
}
\details{
This function calculates weighted estimates for totals and means using suitable weights depending on the class of \code{deskott}: calibrated weights for class \code{kott.cal.design} and direct weights otherwise. Standard errors are calculated using the extended DAGJK method [Kott 99-01].

The mandatory argument \code{y} identifies the variables of interest, that is the variables for which estimates are to be calculated. The corresponding formula must be of the type \code{y=~var1+\ldots+varn}. The \code{deskott} variables referenced by \code{y} must be \code{numeric} or \code{factor} and must not contain any missing value (\code{NA}). It is admissible to specify for \code{y} "mixed" formulas that simultaneously contain quantitative (\code{numeric}) variables and qualitative (\code{factor}) variables.

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{kottby} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. The \code{deskott} variables referenced by \code{by} (if any) must be \code{factor} and must not contain any missing value (\code{NA}). 

The optional argument \code{estimator} makes it possible to select the desired estimator. If\cr \code{estimator="total"} (the default option), \code{kottby} calculates, for a given variable of interest \code{vark}, the estimate of the total (when \code{vark} is \code{numeric}) or the estimate of the absolute frequency distribution (when \code{vark} is \code{factor}). Similarly, if \code{estimator="mean"}, the function calculates the estimate of the mean (when \code{vark} is \code{numeric}) or the the estimate of the relative frequency distribution (when \code{vark} is \code{factor}).

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}.
}
\value{
The return value depends on the value of the input parameters. In the most general case, the function returns an object of class \code{list} (typically a list made up of data frames).
}
\note{
The advantage of the DAGJK method over the traditional jackknife is that, unlike the latter, it remains computationally manageable even when dealing with "complex and big" surveys (tens of thousands of PSUs arranged in a large number of strata with widely varying sizes). In fact, the DAGJK method is known to provide, for a broad range of sampling designs and estimators, (near) unbiased standard error estimates even with a "small" number (e.g. a few tens) of replicate weights. On the other hand, if the number of replicates is not large, it seems defensible to use a \emph{t distribution} (rather than a normal distribution) for calculating the confidence intervals. In line with what was proposed in [Kott 99-01], given an input \code{kott.design} object with \code{nrg} random groups, \code{kottby} builds the confidence intervals making use of a t distribution with \code{nrg-1} degrees of freedom.
}
\author{Diego Zardetto}
\references{
Kott, Phillip S. (1999) \emph{"The Extended Delete-A-Group Jackknife"}. Bulletin of the International Statistical Instititute. 52nd Session. Contributed Papers. Book 2, pp. 167-168.

Kott, Phillip S. (2001) \emph{"The Delete-A-Group Jackknife"}. Journal of Official Statistics, Vol.17, No.4, pp. 521-526.
}
\seealso{
\code{\link{kott.ratio}} for estimating ratios between totals, \code{\link{kott.quantile}} for estimating quantiles, \code{\link{kott.regcoef}} for estimating regression coefficients and \code{\link{kottby.user}} for calculating estimates based on user-defined estimators.
}
\examples{
data(data.examples)

# Creation of a kott.design object:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)


# Estimate of the total of 3 quantitative variables for the whole
# population:
kottby(kdes,~y1+y2+y3)


# Estimate of the total of the same 3 variables by sex: 
kottby(kdes,~y1+y2+y3,~sex)


# Estimate of the mean of the same 3 variables by marstat and sex:
kottby(kdes,~y1+y2+y3,~marstat:sex,estimator="mean")


# Estimate of the absolute frequency distribution of the qualitative
# variable age5c for the whole population:
kottby(kdes,~age5c)


# Estimate of the relative frequency distribution of the qualitative
# variable marstat by sex:
kottby(kdes,~marstat,~sex,estimator="mean")


# The same with confidence intervals at a confidence level of 0.9:
kottby(kdes,~marstat,~sex,estimator="mean",conf.int=TRUE,conf.lev=0.9)


# Quantitative and qualitative variables together: estimate of the
# total for y3 and of the absolute frequency distribution of marstat,
# by sex:
kottby(kdes,~y3+marstat,~sex)


# Lonely PSUs do not give rise to NaNs in the standard errors:
kdes.lpsu<-kottdesign(data=example,ids=~towcod+famcod,strata=~stratum,
           weights=~weight,nrg=15)
kottby(kdes.lpsu,~x1+x2+x3)
}
\keyword{survey}