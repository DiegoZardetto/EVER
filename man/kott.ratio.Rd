\name{kott.ratio}
\alias{kott.ratio}
\title{Estimation of ratios between totals}
\description{
Calculates estimates, standard errors and confidence intervals for ratios between totals in subpopulations.
}
\usage{
kott.ratio(deskott, num, den, by = NULL, 
           vartype = c("se", "cv", "cvpct", "var"), 
           conf.int = FALSE, conf.lev = 0.95)
}
\arguments{
  \item{deskott}{Object of class \code{kott.design} containing the replicated survey data.}
  \item{num}{Formula defining the numerator variables for the ratio estimator.}
  \item{den}{Formula defining the denominator variables for the ratio estimator.}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL} (the default option) estimates refer to the whole population.}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose one or more of: standard error (the default), coefficient of variation, percent coefficient of variation, or variance.}
  \item{conf.int}{Boolean (\code{logical}) value to request confidence intervals for the estimates: the default is \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
}
\details{
This function calculates weighted estimates for ratios between totals of quantitative variables using suitable weights depending on the class of \code{deskott}: calibrated weights for class \code{kott.cal.design} and direct weights otherwise. Standard errors are calculated using the extended DAGJK method [Kott 99-01].

The mandatory argument \code{num} (\code{den}) identifies the variables whose totals appear as the numerator (denominator) in the ratio estimator: the corresponding formula must be of the type \cr
\code{num=~num1+\ldots+numk} (\code{den=~den1+\ldots+denl}). The function calculates estimates for ratios between homologous variables in \code{num} and \code{den}; if \code{num} and \code{den} contain a different number of variables the shortest argument will be tacitly recycled. The \code{deskott} variables referenced by \code{num} (\code{den}) must be \code{numeric} and must not contain any missing value (\code{NA}).

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{kottby} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. The \code{deskott} variables referenced by \code{by} (if any) must be \code{factor} and must not contain any missing value (\code{NA}).

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}. Given an input \code{kott.design} object with \code{nrg} random groups, \code{kott.ratio} builds the confidence intervals making use of a t distribution with \code{nrg-1} degrees of freedom.
}
\value{
The return value depends on the value of the input parameters. In the most general case, the function returns an object of class \code{list} (typically a list made up of data frames).
}
\section{Warning}{
It is possible that, in certain subpopulations, the estimate of the total of some \code{den} variables turns out to be zero for the original sample in \code{deskott} and/or for some of its \code{nrg} replicates. In these cases, \code{kott.ratio} (i) returns \code{NaN} for the estimates and/or for the corresponding standard errors and (ii) prints a \code{warning} message.
}
\author{Diego Zardetto}
\references{
Kott, Phillip S. (1999) \emph{"The Extended Delete-A-Group Jackknife"}. Bulletin of the International Statistical Instititute. 52nd Session. Contributed Papers. Book 2, pp. 167-168.

Kott, Phillip S. (2001) \emph{"The Delete-A-Group Jackknife"}. Journal of Official Statistics, Vol.17, No.4, pp. 521-526.
}
\seealso{
\code{\link{kottby}} for estimating totals and means, \code{\link{kott.quantile}} for estimating quantiles, \code{\link{kott.regcoef}} for estimating regression coefficients and \code{\link{kottby.user}} for calculating estimates based on user-defined estimators.
}
\examples{
data(data.examples)

# Creation of a kott.design object:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)

# Estimate of the ratios y1/x1, y2/x2 e y3/x3 by marstat:
kott.ratio(kdes,~y1+y2+y3,~x1+x2+x3,by=~marstat)

# Estimate of the ratios z/x1, z/x2 e z/x3
# for the whole population (notice the recycling rule):
kott.ratio(kdes,~z,~x1+x2+x3,conf.int=TRUE)

# Estimators of means can be thought as 
# ratio estimators:
kottby(kdes,~income,estimator="mean")  
kott.ratio(kott.addvars(kdes,population=1),num=~income,den=~population)
}
\keyword{survey}