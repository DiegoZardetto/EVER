\name{kott.regcoef}
\alias{kott.regcoef}
\title{Estimation of linear regression coefficients}
\description{
Calculates estimates, standard errors and confidence intervals for regression coefficients in subpopulations.
}
\usage{
kott.regcoef(deskott, model, by = NULL, 
             vartype = c("se", "cv", "cvpct", "var"),
             conf.int = FALSE, conf.lev = 0.95)
}
\arguments{
  \item{deskott}{Object of class \code{kott.design} containing the replicated survey data.}
  \item{model}{Formula giving a symbolic description of the linear model.}
  \item{by}{Formula specifying the variables that define the "estimation domains". If \code{NULL} (the default option) estimates refer to the whole population.}
  \item{vartype}{\code{character} vector specifying the desired variability estimators. It is possible to choose one or more of: standard error (the default), coefficient of variation, percent coefficient of variation, or variance.}
  \item{conf.int}{Boolean (\code{logical}) value to request confidence intervals for the estimates: the default is \code{FALSE}.}
  \item{conf.lev}{Probability specifying the desired confidence level: the default value is \code{0.95}.}
}
\details{
This function calculates weighted estimates of linear regression coefficients using suitable weights depending on the class of \code{deskott}: calibrated weights for class \code{kott.cal.design} and direct weights otherwise. Standard errors are calculated using the extended DAGJK method [Kott 99-01].

The mandatory argument \code{model} specifies, by means of a symbolic \code{\link{formula}}, the linear regression model whose coefficients are to be estimated. \code{model} must have the form \code{response ~ terms} where \code{response} is the (numeric) response variable and \code{terms} represents a series of terms which specifies a linear predictor for \code{response}. Variables referenced by \code{model} must not contain any missing value (\code{NA}).

The optional argument \code{by} specifies the variables that define the "estimation domains", that is the subpopulations for which the estimates are to be calculated. If \code{by=NULL} (the default option), the estimates produced by \code{kottby} refer to the whole population. Estimation domains must be defined by a formula: for example the statement \code{by=~B1:B2} selects as estimation domains the subpopulations determined by crossing the modalities of variables \code{B1} and \code{B2}. The \code{deskott} variables referenced by \code{by} (if any) must be \code{factor} and must not contain any missing value (\code{NA}).

The \code{conf.int} argument allows to request the confidence intervals for the estimates. By default \code{conf.int=FALSE}, that is the confidence intervals are not provided. 

Whenever confidence intervals are requested (i.e. \code{conf.int=TRUE}), the desired confidence level can be specified by means of the \code{conf.lev} argument. The \code{conf.lev} value must represent a probability (\code{0<=conf.lev<=1}) and its default is chosen to be \code{0.95}. Given an input \code{kott.design} object with \code{nrg} random groups and a regression \code{model} with \emph{p} predictors plus an intercept term, \code{kott.regcoef} builds the confidence intervals making use of a t distribution with \code{nrg-p-1} degrees of freedom.
}
\value{
The return value depends on the value of the input parameters. In the most general case, the function returns an object of class \code{list} (typically a list made up of data frames).
}
\author{Diego Zardetto}
\references{
Kott, Phillip S. (1999) \emph{"The Extended Delete-A-Group Jackknife"}. Bulletin of the International Statistical Instititute. 52nd Session. Contributed Papers. Book 2, pp. 167-168.

Kott, Phillip S. (2001) \emph{"The Delete-A-Group Jackknife"}. Journal of Official Statistics, Vol.17, No.4, pp. 521-526.
}
\seealso{
\code{\link{kottby}} for estimating totals and means, \code{\link{kott.ratio}} for estimating ratios between totals, \code{\link{kott.quantile}} for estimating quantiles and \code{\link{kottby.user}} for calculating estimates based on user-defined estimators.
}
\examples{
data(data.examples)

# Creation of a kott.design object:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)


# A model with one predictor and no intercept:
kott.regcoef(kdes,income~z-1)

# ...compare with ratio estimator:
kott.ratio(kott.addvars(kdes,income.mult.z=income*z,z2=z^2),~income.mult.z,~z2)


# A model with a factor term and no intercept:
kott.regcoef(kdes,income~age5c-1)

# ...compare with mean estimator in subpopulations:
kottby(kdes,~income,~age5c,estimator="mean")

# ...and with regression coefficients (for a different model)
# in subpopulations:
kott.regcoef(kdes,income~1,~age5c)


# An awkward model with many coefficients:
kott.regcoef(kdes,income~z:age5c+x3+marstat-1)
}
\keyword{survey}