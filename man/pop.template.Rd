\name{pop.template}
\alias{pop.template}
\title{Template data frame for known population totals}
\description{
Construct a \emph{"template"} data frame to store known population totals for a calibration problem.
}
\usage{
pop.template(data, calmodel, partition = FALSE)
}
\arguments{
  \item{data}{Data frame of survey data.}
  \item{calmodel}{Formula defining the linear structure of the calibration model.}
  \item{partition}{Formula specifying the variables that define the "calibration domains" for the model. \code{FALSE} (the default) implies no calibration domains.}
}
\details{
This function creates an object of class \code{pop.totals}. A \code{pop.totals} object is made up by the union of a data frame whose structure conforms to the standard required by \code{kottcalibrate} for the known totals and the metadata describing the calibration problem.

The mandatory argument \code{data} must identify the survey data frame on which the calibration problem is defined (or, as an alternative, a \code{kott.design} object built upon that data frame). 

The mandatory argument \code{calmodel} symbolically defines the calibration model you want to use: it identifies the auxiliary variables and the constraints for the calibration problem. The \code{deskott} variables referenced by \code{calmodel} must be \code{numeric} or \code{factor} and must not contain any missing value (\code{NA}).

The optional argument \code{partition} specifies the variables that define the calibration domains for the model. The default value (\code{FALSE}) means either that there are not calibration domains or that you want to solve the problem globally (even though it could be factorised). If a formula is passed through the \code{partition} argument the program checks that \code{calmodel} actually describes a "reduced model", that is it does not reference any of the partition variables; if this is not the case, the program stops and prints an error message. The \code{deskott} variables referenced by \code{partition} (if any) must be \code{factor} and must not contain any missing value (\code{NA}).
}
\value{
An object of class \code{pop.totals}. The data frame it contains is a \emph{"template"} in the sense that all the known totals it must be able to store are missing (\code{NA}). However, this data frame has a structure that complies with the standard required by \code{kottcalibrate} (provided the latter is invoked with the same \code{calmodel} and \code{partition} values used to create the template).\cr
The operation of filling the template's \code{NA}s with the actual values of the corresponding population totals is, obviously, to be done by the user.

The \code{pop.totals} class is a specialisation of the \code{data.frame} class; this means that an object built by \code{pop.template} inherits from the \code{data.frame} class and you can use on it every method defined on that class.
}
\author{Diego Zardetto}
\seealso{
\code{\link{kottcalibrate}} for calibrating replicate weights, \code{\link{population.check}} to check that the known totals data frame satisfies the standard required by \code{kottcalibrate}.
}
\examples{
# Creation of population totals template data frames for different
# calibration problems (if the calibration models can be factorised 
# both a global and an iterative solution are given):

data(data.examples)

# 1) Calibration on the total number of units in the population:
pop.template(data=example,calmodel=~1)


# 2) Calibration on the total number of units in the population
#    and on the marginal distribution of marstat (notice that the
#    total for the first level "married" of the marstat factor
#    variable is missing because it can be deduced from
#    the remaining totals):
pop.template(data=example,calmodel=~marstat)


# 3) Calibration on the marginal distribution of marstat (you
#    must explicitly remove the intercept term in the 
#    calibration model adding -1 to the calmodel formula):
pop.template(data=example,calmodel=~marstat-1)


# 4) Calibration (global solution) on the joint distribution of sex
#    and marstat:
pop.template(data=example,calmodel=~sex:marstat-1)

# 4.1) Calibration (iterative solution) on the joint distribution
#      of sex and marstat:
#      4.1.1) Using sex to define calibration domains:
pop.template(data=example,calmodel=~marstat-1,partition=~sex)

#      4.1.2) Using marstat to define calibration domains:
pop.template(data=example,calmodel=~sex-1,partition=~marstat)


# 5) Calibration (global solution) on the total for the quantitative
#    variable x1 and on the marginal distribution of the qualitative
#    variable age5c, in the subpopulations defined by crossing sex
#    and marstat:
pop.template(data=example,calmodel=~(age5c+x1-1):sex:marstat)

# 5.1) The same problem with iterative solutions:
#      5.1.1) Using sex to define calibration domains:
pop.template(data=example,calmodel=~(age5c+x1-1):marstat,partition=~sex)

#      5.1.2) Using marstat to define calibration domains:
pop.template(data=example,calmodel=~(age5c+x1-1):sex,partition=~marstat)

#      5.1.3) Using sex and marstat to define calibration domains:
pop.template(data=example,calmodel=~age5c+x1-1,partition=~sex:marstat)
}
\keyword{survey}