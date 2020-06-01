\name{kottcalibrate}
\alias{kottcalibrate}
\alias{kottcal.status}
\title{Calibration of replicate weights}
\description{
Adds to a \code{kott.design} object the calibrated weights columns (one for each replicate weight, plus one for the initial weights). 
}
\usage{
kottcalibrate(deskott, df.population,
        calmodel = if (inherits(df.population, "pop.totals"))
                       attr(df.population, "calmodel"),
        partition = if (inherits(df.population, "pop.totals"))
                        attr(df.population, "partition") else FALSE,
        calfun = c("linear", "raking", "logit"),
        bounds = c(-Inf, Inf), aggregate.stage = NULL, maxit = 50,
        epsilon = 1e-07, force.rep = FALSE)
}
\arguments{
  \item{deskott}{Object of class \code{kott.design} containing the replicated survey data.}
  \item{df.population}{Data frame containing the known population totals for the auxiliary variables.}
  \item{calmodel}{Formula defining the linear structure of the calibration model.}
  \item{partition}{Formula specifying the variables that define the "calibration domains" for the model (see 'Details'); \code{FALSE} (the default) implies no calibration domains.}
  \item{calfun}{\code{character} specifying the distance function for the calibration process; the default is \code{"linear"}.}  
  \item{bounds}{Allowed range for the ratios between calibrated and initial weights; the default is \code{c(-Inf,Inf)}.}  
  \item{aggregate.stage}{An integer: if specified, causes the calibrated weights to be constant within sampling units at this stage.}
  \item{maxit}{Maximum number of iterations for the Newton-Raphson algorithm; the default is \code{50}.}
  \item{epsilon}{Tolerance for the relative differences between the population totals and the corresponding estimates based on the claibrated weights; the default is \code{10^-7}.}
  \item{force.rep}{If \code{TRUE}, whenever the calibration algorithm does not converge for a given set of replicate weights, forces the function to return a value (see 'Details'); the default is \code{FALSE}.} 
}
\details{
This function creates an object of class \code{kott.cal.design}. A \code{kott.cal.design} object is made up by the union of the (calibrated) replicated survey data and the metadata describing the sampling design. \code{kott.cal.design} objects make it possible to estimate the variance of calibration estimators [Deville, Sarndal 92] using the extended \emph{"Delete-A-Group Jackknife"} method [Kott 2008].

The mandatory argument \code{calmodel} symbolically defines the calibration model you want to use, that is - in the language of the generalised regression estimator - the assisting linear regression model underlying the calibration problem [Wilkinson, Rogers 73]. More specifically, the \code{calmodel} formula identifies the auxiliary variables and the constraints for the calibration problem. For example, \code{calmodel=~(X+Z):C+(A+B):D} defines the calibration problem in which constraints are imposed: (i) on the auxiliary (quantitative) variables \code{X} and \code{Z} within the subpopulations identified by the (qualitative) classification variable \code{C} and, at the same time, (ii) on the absolute frequency of the (qualitative) variables \code{A} and \code{B} within the subpopulations identified by the (qualitative) classification variable \code{D}.\cr
The \code{deskott} variables referenced by \code{calmodel} must be \code{numeric} or \code{factor} and must not contain any missing value (\code{NA}).

Problems for which one or more qualitative variables can be \emph{"factorised"} in the formula that specifies the calibration model, are particularly interesting. These variables split the population into non-overlapping subpopulations known as \emph{"calibration domains"} for the model. An example is provided by the statement \code{calmodel=~(A+B+X+Z):D} in which the variable that identifies the calibration domains is \code{D}; similarly, the formula \code{calmodel=~(A+B+X+Z):D1:D2} identifies as calibration domains the subpopulations determined by crossing the modalities of \code{D1} and \code{D2}. The interest in models of this kind lies in the fact that the \emph{global} calibration problem they describe can, actually, be broken down into \emph{local} subproblems, one per calibration domain, which can be solved separately [Vanderhoeft 01]. Thus, for example, the global problem defined by \code{calmodel=~(A+B+X+Z):D} is equivalent to the sequence of problems defined by the \emph{"reduced model"} \code{calmodel=~A+B+X+Z} in each of the domains identified by the modalities of \code{D}. The opportunity to separately solve the subproblems related to different calibration domains achieves a significant reduction in computation complexity: the gain increases with increasing survey data size and (most importantly) with increasing auxiliary variables number.

The optional argument \code{partition} makes it possible to choose, in cases in which the calibration problem can be factorised, whether to solve the problem globally or iteratively (that is, separately for each calibration domain). The global solution (which is the default option) can be selected invoking the \code{kottcalibrate} function with \code{partition=FALSE}. To request the iterative solution - a strongly recommended option when dealing with a lot of auxiliary variables and big data sizes - it is necessary to specify via \code{partition} the variables defining the calibration domains for the model. If a formula is passed through the \code{partition} argument (for example: \code{partition=~D1:D2}), the program checks that \code{calmodel} actually describes a "reduced model" (for example: \code{calmodel=~X+Z+A+B}), that is it does not reference any of the partition variables; if this is not the case, the program stops and prints an error message.\cr
The \code{deskott} variables referenced by \code{partition} (if any) must be \code{factor} and must not contain any missing value (\code{NA}).

The mandatory argument \code{df.population} is used to specify the known totals of the auxiliary variables referenced by \code{calmodel} within the subpopulations (if any) identified by \code{partition}. These known totals must be stored in a data frame whose structure (i) depends on the values of \code{calmodel} and \code{partition} and (ii) must conform to a standard. In order to facilitate understanding of and compliance with this standard, the \pkg{EVER} package provides the user with two functions: \code{\link{pop.template}} and \code{\link{population.check}}. The \code{pop.template} function is able to guide the user in constructing the known totals data frame for a specific calibration problem, while the \code{population.check} function allows to check whether a known totals data frame conforms to the standard required by \code{kottcalibrate}. In any case, if the \code{df.population} data frame does not comply with the standard, the \code{kottcalibrate} function stops and prints an error message: the meaning of the message should help the user diagnose the cause of the problem.

The \code{calfun} argument identifies the distance function to be used in the calibration process. Three built-in functions are provided: \code{"linear"}, \code{"raking"}, and \code{"logit"}. The default is \code{"linear"}, which corresponds to the euclidean metric.

The \code{bounds} argument allows to add \emph{"range constraints"} to the calibration problem. To be precise, the interval defined by \code{bounds} will contain the values of the ratios between final (calibrated) and initial (direct) weights. The default value is \code{c(-Inf,Inf)}, i.e. no range constraints are imposed. These constraints are optional unless the \code{"logit"} function is selected: in the latter case the range defined by \code{bounds} has to be finite.       

The value passed by the \code{aggregate.stage} argument must be an integer between 1 and the number of sampling stages of \code{deskott}. If specified, causes the calibrated weights to be constant within sampling units selected at the \code{aggregate.stage} stage (actually this is only ensured if the initial weights had already this property, as is sometimes the case in multistage cluster sampling). If not specified, the calibrated weights may differ even for sampling units with identical initial weights. The same holds if some final units belonging to the same cluster selected at the stage \code{aggregate.stage} fall in distinct calibration domains (i.e. if the domains defined by \code{partition} "cut across" the \code{aggregate.stage}-stage clusters).

The \code{maxit} argument sets the maximum number of iteration for the Newton-Raphson algorithm that is used to solve the calibration problem. The default value is \code{50}.

The \code{epsilon} argument determines the convergence criterion for the optimisation algorithm: it fixes the maximum allowed value for the relative differences between the population totals and the corresponding estimates based on the claibrated weights. The default value is \code{10^-7}.

If the number of replicates for \code{deskott} (the input object of class \code{kott.design}) is \code{nrg}, the function \code{kottcalibrate} is in charge of solving \code{nrg+1} distinct calibration problems. In fact, the calibrated weights calculated by \code{kottcalibrate} must ensure that the known population totals are exactly reproduced not only by the original sample, but also by all its \code{nrg} replicates. Should this requirement fail, the DAGJK method would end up with a biased variance estimator [Kott 2008]. It is, however, possible (more likely when range constraints are imposed) that, for some of the \code{nrg+1} distinct calibration problems and for the given values of \code{epsilon} and \code{maxit}, the solving algorithm does not converge. In this case \code{kottcalibrate} by default stops and prints an error message. On the contrary if \code{force.rep = TRUE}, \emph{provided that the failure to converge pertains only to the replicate weights}, the function is forced to return the best approximation achieved for the corresponding calibrated weights. When this occurs, DAGJK standard errors estimates built on the object returned by \code{kottcalibrate} will be biased.
}
\section{Calibration process diagnostics}{
If the number of replicates for \code{deskott} is \code{nrg}, the function \code{kottcalibrate} is in charge of solving \code{nrg+1} distinct calibration problems. When, dealing with a factorisable calibration problem, the user selects the iterative solution, each one of the above mentioned problems is split into as many \emph{sub-problems} as the number of subpopulations defined by \code{partition}. A calibration process with such a complex structure needs some ad hoc tool for error diagnostics. For this purpose, every call to \code{kottcalibrate} creates, by side effect, a dedicated data structure named \code{kottcal.status} into the \code{.GlobalEnv}. \code{kottcal.status} is a list with two components: the first, \code{"call"}, identifies the call to \code{kottcalibrate} that generated the list, the second, \code{return.code}, is a matrix each element of which identifies the return code of a specific calibration sub-problem. The meaning of the return codes is as follows:
\describe{
\item{\code{-1}}{not yet tackled sub-problem;}
\item{\code{0}}{solved sub-problem (convergence achieved);}
\item{\code{1}}{unsolved sub-problem (no convergence): output forced.}
}
Recall that the latter return code may only occur if \code{force.rep = TRUE}.\cr
In case of error, users can exploit \code{kottcal.status} to identify the sub-problem from which the error stemmed, hence taking a step forward to eliminate it.
}
\value{
An object of class \code{kott.cal.design}. The data frame it contains includes (in addition to the data already stored in \code{deskott}) the calibrated weights columns (one for each replicate weight, plus one for the initial weights, \code{nrg+1} in all). The names of these columns are obtained by pasting the name of the initial weights column with the string \code{".cal"} and the indices \code{NULL, 1, 2, \ldots, nrg}.\cr
The \code{kott.cal.design} class is a specialisation of the \code{kott.design} class; this means that an object created by \code{kottcalibrate} inherits from the \code{data.frame} class and you can use on it every method defined on that class.
}
\author{Diego Zardetto}
\references{
Deville, J.C., Sarndal, C.E. (1992) \emph{"Calibration Estimators in Survey Sampling"}, Journal of the American Statistical Association, Vol. 87, No. 418, pp.376-382.

Kott, Phillip S. (2008) \emph{"Building a Better Delete-a-Group Jackknife for a Calibration Estimator"}, NASS Research Report, NASS: Washington, DC.

Wilkinson, G.N., Rogers, C.E. (1973) \emph{"Symbolic Description of Factorial Models for Analysis of Variance"}, Journal of the Royal Statistical Society, series C (Applied Statistics), Vol. 22, pp. 181-191.

Vanderhoeft, C. (2001) \emph{"Generalized Calibration at Statistic Belgium"}, Statistics Belgium Working Paper n. 3, \url{http://www.statbel.fgov.be/studies/paper03_en.asp}.

Lumley, T. (2006) \emph{"survey: analysis of complex survey samples"}, \url{http://cran.at.r-project.org/web/packages/survey/index.html}.

Scannapieco, M., Zardetto, D., Barcaroli, G. (2007) \emph{"La Calibrazione dei Dati con R: una Sperimentazione sull'Indagine Forze di Lavoro ed un Confronto con GENESEES/SAS"}, Contributi Istat n. 4., \url{http://www.istat.it/dati/pubbsci/contributi/Contributi/contr_2007/2007_4.pdf}.
}
\seealso{
\code{\link{desc}} for a concise description of \code{kott.design} objects, \code{\link{kottby}}, \code{\link{kott.ratio}}, \code{\link{kott.regcoef}}, \code{\link{kott.quantile}} and \code{\link{kottby.user}} for calculating estimates and standard errors, \code{\link{pop.template}} for constructing known totals data frames in compliance with the standard required by \code{kottcalibrate}, \code{\link{population.check}} to check that the known totals data frame satisfies that standard, \code{\link{bounds.hint}} to obtain an hint for range restricted calibration. 
}
\examples{
# Calibration of a kott.design object according to different calibration
# models (the known totals data frames pop01, \ldots, pop05p and the bounds
# vector are contained in the data.examples file).
# For the examples relating to calibration models that can be factorised 
# both a global and an iterative solution are given.

data(data.examples)

# Creation of the object to be calibrated:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)


# 1) Calibration on the total number of units in the population
#   (totals in pop01):
kdescal01<-kottcalibrate(deskott=kdes,df.population=pop01,calmodel=~1,
           calfun="logit",bounds=bounds,aggregate.stage=2)

# Checking the result (the function 'ones' is contained
# in the data.examples file):
kottby.user(kdescal01,user.estimator=ones)


# 2) Calibration on the marginal distributions of sex and marstat
#    (totals in pop02):
kdescal02<-kottcalibrate(deskott=kdes,df.population=pop02,
           calmodel=~sex+marstat-1,calfun="logit",bounds=bounds,
           aggregate.stage=2)

# Checking the result:
kottby(kdescal02,~sex+marstat)


# 3) Calibration (global solution) on the joint distribution of sex
#    and marstat (totals in pop03):
kdescal03<-kottcalibrate(deskott=kdes,df.population=pop03,
           calmodel=~marstat:sex-1,calfun="logit",bounds=bounds)

# Checking the result:
kottby(kdescal03,~sex,~marstat) # or: kottby(kdescal03,~marstat,~sex)

# which, obviously, is not respected by kdescal02 (notice the size of SE):
kottby(kdescal02,~sex,~marstat)


# 3.1) Again a calibration on the joint distribution of sex and marstat
#      but, this time, with the iterative solution (partition=~sex,
#      totals in pop03p):
kdescal03p<-kottcalibrate(deskott=kdes,df.population=pop03p,
            calmodel=~marstat-1,partition=~sex,calfun="logit",
            bounds=bounds)

# Checking the result:
kottby(kdescal03p,~sex,~marstat)


# 4) Calibration (global solution) on the totals for the quantitative
#    variables x1, x2 and x3 in the subpopulations defined by the
#    regcod variable (totals in pop04):
kdescal04<-kottcalibrate(deskott=kdes,df.population=pop04,
           calmodel=~(x1+x2+x3-1):regcod,calfun="logit",
           bounds=bounds,aggregate.stage=2)

# Checking the result:
kottby(kdescal04,~x1+x2+x3,~regcod)


# 4.1) Same problem with the iterative solution (partition=~regcod,
#      totals in pop04p):
kdescal04p<-kottcalibrate(deskott=kdes,df.population=pop04p,
            calmodel=~x1+x2+x3-1,partition=~regcod,calfun="logit",
            bounds=bounds,aggregate.stage=2)

# Checking the result:
kottby(kdescal04p,~x1+x2+x3,~regcod)


# 5) Calibration (global solution) on the total for the quantitative
#    variable x1 and on the marginal distribution of the qualitative
#    variable age5c, in the subpopulations defined by crossing sex
#    and marstat (totals in pop05):
kdescal05<-kottcalibrate(deskott=kdes,df.population=pop05,
           calmodel=~(age5c+x1-1):sex:marstat,calfun="logit",
           bounds=bounds,force.rep=TRUE)

# Calibration process diagnostics:
kottcal.status

# Checking the result:
kottby(kdescal05,~age5c+x1,~sex:marstat)


# 5.1) Same problem with the iterative solution (partition=~sex:marstat,
#      totals in pop05p):
kdescal05p<-kottcalibrate(deskott=kdes,df.population=pop05p,
            calmodel=~age5c+x1-1,partition=~sex:marstat,
            calfun="logit",bounds=bounds,force.rep=TRUE)

# Calibration process diagnostics:
kottcal.status

# Checking the result:
kottby(kdescal05p,~age5c+x1,~sex:marstat)


# Notice that 3.1 e 5.1) do not impose the aggregate.stage=2
# condition. This condition cannot, in fact, be fulfilled because
# in both cases the domains defined by partition "cut across"
# the kdes second stage clusters (households). To compare the results,
# the same choice was also made for 3) e 5).
}
\keyword{survey}