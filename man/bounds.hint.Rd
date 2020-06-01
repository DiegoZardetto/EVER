\name{bounds.hint}
\alias{bounds.hint}
\title{A hint for range restricted calibration}
\description{
Suggests a sound \code{bounds} value for which \code{kottcalibrate} is likely to converge.
}
\usage{
bounds.hint(deskott, df.population, 
    calmodel = if (inherits(df.population, "pop.totals"))
                   attr(df.population, "calmodel"), 
    partition = if (inherits(df.population, "pop.totals")) 
                    attr(df.population, "partition") else FALSE)
}
\arguments{
  \item{deskott}{Object of class \code{kott.design} containing the replicated survey data.}
  \item{df.population}{Data frame containing the known population totals for the auxiliary variables.}
  \item{calmodel}{Formula defining the linear structure of the calibration model.}
  \item{partition}{Formula specifying the variables that define the "calibration domains" for the model; \code{FALSE} (the default) implies no calibration domains.}
}
\details{
The function \code{bounds.hint} returns a \code{bounds} value for which \code{kottcalibtrate} is \emph{likely} to converge. This interval is just a sound hint, \emph{not} an exact result (see 'Note').

The mandatory argument \code{deskott} identifies the \code{kott.design} object on which the calibration problem is defined.

The mandatory argument \code{df.population} identifies the known totals data frame.

The argument \code{calmodel} symbolically defines the calibration model you want to use: it identifies the auxiliary variables and the constraints for the calibration problem. The \code{deskott} variables referenced by \code{calmodel} must be \code{numeric} or \code{factor} and must not contain any missing value (\code{NA}). The argument can be omitted provided \code{df.population} is an object of class \code{pop.totals} (see \code{\link{population.check}}). 

The optional argument \code{partition} specifies the variables that define the calibration domains for the model. The default value (\code{FALSE}) means either that there are not calibration domains or that you want to solve the problem globally (even though it could be factorised). The \code{deskott} variables referenced by \code{partition} (if any) must be \code{factor} and must not contain any missing value (\code{NA}). The argument can be omitted provided \code{df.population} is an object of class \code{pop.totals} (see \code{\link{population.check}}).
}
\value{
A numeric vector of length 2, representing the \emph{suggested} value for the \code{bounds} argument of \code{kottcalibrate}. The attributes of that vector store additional information, which can lead to better understand why a given calibration problem is (un)feasible (see 'Examples').
}
\note{
Assessing the feasibility of an arbitrary calibration problem is not an easy task. The problem is even more difficult whenever additional \emph{"range restrictions"} are imposed. Indeed, even if one assumes that the calibration constraints define a consistent system, one also has to choose the \code{bounds} such that the feasible region is non-empty.

One can argue that there must exist a minimun-length interval \eqn{I=[L,U]} such that, if it is covered by \code{bounds}, the specified calibration problem is feasible. Unfortunately in order to compute exactly that minimun-length interval \eqn{I} one should solve a big linear programming problem [Vanderhoeft 01]. As an alternative, a trial and error procedure has been frequently proposed [Deville et al 1993; Sautory 1993]: (i) start with a very large interval \code{bounds.0}; (ii) if convergence is achieved, shrink it so as to obtain a new inteval \code{bounds.1}; (iii) repeat until you get a sufficiently tight feasible interval \code{bounds.n}. The drawback is that this procedure can cost a lot of computer time since, for each choice of the \code{bounds}, the full calibration problem has to be solved.

A rather easy task is, on the contrary, the one of finding at least a given specific interval \eqn{I^*=[L^*,U^*]} such that, if it is \emph{not} covered by \code{bounds}, the current calibration problem is \emph{surely unfeasible}. This means that any feasible \code{bounds} value must necessarily contain the \eqn{I^*} interval. The function \code{bounds.hint}: (i) first identifies such an \eqn{I^*} interval (by computing the range of the ratios between known population totals and corresponding direct Horvitz-Thompson estimates), (ii) then builds a new interval \eqn{I^{sugg}} with same midpoint and double length. The latter is the \emph{suggested} value for the \code{bounds} argument of \code{kottcalibrate}. The return value of \code{bounds.hint} should be understood as a useful starting guess for \code{bounds}, even though there is definitely no warranty that the calibration algorithm will actually converge.
}
\author{Diego Zardetto}
\references{
Vanderhoeft, C. (2001) \emph{"Generalized Calibration at Statistic Belgium"}, Statistics Belgium Working Paper n. 3, \url{http://www.statbel.fgov.be/studies/paper03_en.asp}.

Deville, J.C., Sarndal, C.E. and Sautory, O. (1993) \emph{"Generalized Raking Procedures in Survey Sampling"}, Journal of the American Statistical Association, Vol. 88, No. 423, pp.1013-1020.

Sautory, O. (1993) \emph{"La macro CALMAR: Redressement d'un Echantillon par Calage sur Marges"}, Document de travail de la Direction des Statistiques Demographiques et Sociales, no. F9310.
}
\seealso{
\code{\link{kottcalibrate}} for calibrating replicate weights, \code{\link{pop.template}} for constructing known totals data frames in compliance with the standard required by \code{kottcalibrate}, \code{\link{population.check}} to check that the known totals data frame satisfies that standard and \code{\link{g.range}} to compute the range of the obtained g-weights. 
}
\examples{
# Load sample data (the only reason for fixing
# the RNG seed is to achieve reproducible examples)
data(data.examples)
set.seed(123)

# Creation of the object to be calibrated:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)

# Calibration (global solution) on the joint distribution
# of sex and marstat (totals in pop03). Get a hint for feasible bounds:
hint<-bounds.hint(kdes,pop03,~marstat:sex-1)

# Let's first verify if calibration converges with the suggested 
# value for the bounds argument (i.e. c(0.909, 1.062) ):
kdescal03<-kottcalibrate(deskott=kdes,df.population=pop03,
           calmodel=~marstat:sex-1,calfun="logit",bounds=hint)

# Now let's verify that calibration fails, if bounds don't cover
# the interval [0.947, 1.023]:
\dontrun{
kdescal03<-kottcalibrate(deskott=kdes,df.population=pop03,
           calmodel=~marstat:sex-1,calfun="logit",bounds=c(0.95, 1.03))
}


# Calibration (iterative solution) on the totals for the quantitative
# variables x1, x2 and x3 in the subpopulations defined by the
# regcod variable (totals in pop04p): Get a hint for feasible bounds:
hint<-bounds.hint(kdes,pop04p,~x1+x2+x3-1,~regcod)

# Let's verify if calibration converges with the suggested 
# value for the bounds argument (i.e. c(0.133, 2.497) ):
kdescal04p<-kottcalibrate(deskott=kdes,df.population=pop04p,
            calmodel=~x1+x2+x3-1,partition=~regcod,calfun="logit",
            bounds=hint,aggregate.stage=2)
			
# Now let's verify that calibration fails, if bounds don't cover
# the interval [0.724, 1.906]:
\dontrun{
kdescal04p<-kottcalibrate(deskott=kdes,df.population=pop04p,
            calmodel=~x1+x2+x3-1,partition=~regcod,calfun="logit",
            bounds=c(0.71,1.89),aggregate.stage=2)
}

# By analysing kottcal.status one understands that calibration
# failed due to the sub-task identified by replicate.12 and 
# regcod 6:
kottcal.status

# this is easily explained by inspectioning the "bounds"
# attribute of the bounds.hint output object:
hint

# indeed the specified upper bound (1.89) was too low
# for replicate.12 and regcod 6
}
\keyword{survey}