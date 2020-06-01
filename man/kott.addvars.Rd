\name{kott.addvars}
\alias{kott.addvars}
\title{Add variables to replicated objects}
\description{
Modifies a kott.design object by adding new variables to it. 
}
\usage{
kott.addvars(deskott, ...)
}
\arguments{
  \item{deskott}{Object of class \code{kott.design} containing the replicated survey data.}
  \item{\ldots}{\code{tag = expr} arguments defining columns to be added to \code{deskott}.}
}
\details{
This function adds to the data frame contained in \code{deskott} the \emph{new} variables defined by the \code{tag = expr} arguments. A \code{tag} can be specified either by means of an identifier or by a character string; \code{expr} can be any expression that it makes sense to evaluate in the \code{deskott} environment.

For each argument \code{tag = expr} bound to the formal argument \code{\ldots} the added column will have \emph{name} given by the \code{tag} value and \emph{values} obtained by evaluating the \code{expr} expression on \code{deskott}. Any input expression unsupplied with a \code{tag} will be ignored and will therefore have no effect on the \code{kott.addvars} return value.

Variables to be added to the input replicated object have to be \emph{new}: namely it is not possible to use \code{kott.addvars} to modify the values in a pre-existing \code{deskott} column.
}
\value{
An object of the same class of \code{deskott}, containing new variables but supplied with exactly the same metadata.
}
\author{Diego Zardetto}
\examples{
data(data.examples)

# Creation of a kott.design object:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)

# Adding the new 'population' variable to estimate the number
# of final units in the population: 
kdes2<-kott.addvars(kdes,population=1)
kottby(kdes2,~population)

# Recoding a qualitative variable:
kdes2<-kott.addvars(kdes,agerange=as.factor(ifelse(age5c==1,
                                            "young","not-young")))
kottby(kdes2,~agerange,est="mean")
kottby(kdes2,~income,~agerange,estimator="mean",conf.int=TRUE)

# Algebraic operations on numeric variables:
kdes2<-kott.addvars(kdes,q=income/z^2)
kottby(kdes2,~q)
}
\keyword{survey}