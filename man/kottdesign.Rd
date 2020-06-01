\name{kottdesign}
\alias{kottdesign}
\title{Delete-A-Group Jackknife replication}
\description{
Adds to a data frame of survey data the replicate weights calculated according to the \emph{"Delete-A-Group Jackknife"} (DAGJK) method.}
\usage{
kottdesign(data, ids, strata = FALSE, weights, nrg, 
           self.rep.str = FALSE, check.data = FALSE, 
           aux = FALSE)
}
\arguments{
  \item{data}{Data frame of survey data.}
  \item{ids}{Formula identifying clusters selected at subsequent sampling stages (PSUs, SSUs, \ldots).}
  \item{strata}{Formula identifying the stratification variable; \code{FALSE} (the default) implies no stratification.}  
  \item{weights}{Formula identifying the initial weights for the sampling units.}
  \item{nrg}{Number of "random groups" (and replicate weights) you want to create.} 
  \item{self.rep.str}{Formula identifying self-representing strata (SR), if any; \code{FALSE} (the default)
  means no SR strata.}
  \item{check.data}{Boolean (\code{logical}) value to check the correct nesting of \code{data} clusters; the default is \code{FALSE}.}
  \item{aux}{If \code{TRUE} adds columns of auxiliary information to the output data frame.}
}
\details{
This function creates an object of class \code{kott.design}. A \code{kott.design} object is made up by the union of the replicated survey data and the metadata describing the sampling design. The metadata (stored as attributes of the object) are used to enable and guide processing and analyses provided by other functions in the \pkg{EVER} package (such as \code{kottcalibrate}, \code{kottby}, \code{desc}, \ldots).

The \code{data}, \code{ids}, \code{weights} and \code{nrg} arguments are mandatory, while \code{strata}, \code{check.data} and \code{aux} arguments are optional. The \code{data} variables that are referenced by \code{ids}, \code{weights} and \code{strata} (if specified) must not contain any missing value (\code{NA}). 

The \code{ids} argument specifies the cluster identifiers. It is possible to specify a multi-stage sampling design by simply using a formula with the identifiers of clusters selected at subsequent sampling stages. For example, \code{ids=~id.PSU+id.SSU} declares a two-stage sampling in which the first stage units are identified by the \code{id.PSU} variable and second stage ones by the \code{id.SSU} variable. 

The \code{strata} argument identifies the stratification variable. The \code{data} variable referenced by \code{strata} (if specified) must be a \code{factor}. By default the sample is assumed to be non-stratified.

The \code{weights} argument identifies the initial (or direct) weights for the units included in the sample. The \code{data} variable referenced by \code{weights} must be \code{numeric}.

The \code{nrg} argument selects the number of "random groups" (and replicate weights) you want to create by means of the DAGJK method [Kott 98-99-01]. The value of \code{nrg} must be greater than 1 and less than or equal to the number of sampled PSUs (otherwise the function stops and prints an error message). If \code{nrg} equals the number of sampled PSUs, the DAGJK method "reduces" to (that is, it provides identical results to) the traditional stratified jackknife method. The advantage of the DAGJK method over the traditional jackknife is that, unlike the latter, it remains computationally manageable even when dealing with "complex and big" surveys (tens of thousands of PSUs arranged in a large number of strata with widely varying sizes). In fact, the DAGJK method is known to provide, for a broad range of sampling designs and estimators, (near) unbiased standard error estimates even with a "small" number (e.g. a few tens) of replicate weights.

When dealing with a multistage, stratified sampling design that includes \emph{self-representing (SR) strata} (i.e. strata containing PSUs selected with probability 1), the main contribution to the variance of the SR strata arises from the second stage units (\emph{"variance PSUs"}). In this instance, the user can exploit the \code{self.rep.str} argument to specify, by a formula, the \code{data} variable identifying the SR strata: as a result the function will build the variance PSUs and take care of them. When choosing this option, the user must ensure that the variable referenced by \code{self.rep.str} is \code{logical} (with value \code{TRUE} for SR strata and \code{FALSE} otherwise) or \code{numeric} (with value \code{1} for SR strata and \code{0} otherwise).\cr
As an alternative, the user can attend to develop by himself the appropriate identifiers for the sampling units in \code{ids}. To be precise, the identifier for the PSUs (say \code{id.PSU}) must have, in the SR strata, values in correspondence 1:1 (for example they can be equal, provided this does not cause undesired duplications) with those of the SSUs identifier (say \code{id.SSU}).

The optional argument \code{check.data} allows to check the correct nesting of \code{data} clusters (PSUs, SSUs, \ldots). If \code{check.data=TRUE} the function checks that every unit selected at stage \code{k+1} is associated to one and only one unit selected at stage \code{k}. For a stratified design the function checks also the correct nesting of clusters within strata.

The optional argument \code{aux} can usually be ignored: its default value selects the standard behaviour of the function. Invoking \code{kottdesign} with \code{aux=TRUE} can, on the other hand, prove useful for any user who wants to fully understand how the DAGJK method builds the replicate weights. If \code{aux=TRUE}, the output data frame contains auxiliary columns that provide: the number of PSUs per stratum, the number of PSUs per stratum and random group and the multiplicative coefficients that transform the initial weights into replicate weights.
}
\value{
An object of class \code{kott.design}. The data frame it contains includes (in addition to the original survey data):
\item{-}{A new column named \emph{rgi} (random group index) giving the random group to which each sample unit belongs.}
\item{-}{The replicate weights columns (one per random group, \code{nrg} in all), the names of which are obtained by pasting the name of the initial weights column with the indices 1, 2, \ldots, \code{nrg}.}
The \code{kott.design} class is a specialisation of the \code{data.frame} class; this means that an object created by \code{kottdesign} inherits from the \code{data.frame} class and you can use on it every method defined on that class.
}
\note{
The \pkg{EVER} package implements the extended version of the DAGJK method [Kott 99-01]. It guarantees unbiased estimates of standard errors even when the number of PSUs sampled in some strata is small (that is, less than \code{nrg}).

The rigorous [Kott 98-99-01] results were derived under the hypothesis of with replacement selection of PSUs. This means that the DAGJK method cannot include finite population corrections (\emph{fpc}): this restriction is fully reflected in the \pkg{EVER} package.

If only one PSU (\emph{lonely PSU}) has been selected in some non-self-representative strata (NSR), the \code{kottdesign} function does not report an error message, rather a warning one. In fact, the extended DAGJK method automatically removes the contribution of strata containing lonely PSUs from the estimation of standard errors (obviously, their contribution remains when calculating the estimates). This is all the users have to remember, if they come across the warning message produced by \code{kottdesign}. Whenever the described behaviour seems to be undesirable, a viable alternative in order to eliminate the lonely PSUs is to collapse strata in a suitable manner. In such a case, the price to pay is the possibility of ending up with an over-estimation of the standard errors. As far as the strata collapsing strategie is concerned, the \pkg{EVER} package does not provide (in the current version) any support to the user.

Unlike the conventional jackknife method, the DAGJK is a stochastic replication method. If, having fixed the sampling design and the number of replicates, it is applied a number of times to the same sample data frame, generally a different random groups composition results. This means that repeated invocations of the \code{kottdesign} function, even if run with identical actual parameters, generate different \code{kott.design} objects (and, consequently, different standard error estimates). What has been stated obviously does not apply when \code{nrg} equals the number of sampled PSUs. If you really need it, you can however generate exactly the same results for subsequent applications of \code{kottdesign}: you have only to keep fixed the seed of \R's random numbers generator (using the \code{\link{set.seed}} function).
}
\author{Diego Zardetto.}
\references{
Kott, Phillip S. (1998) \emph{"Using the Delete-A-Group Jackknife Variance Estimator in NASS Surveys"}, RD Research Report No. RD-98-01, USDA, NASS: Washington, DC.

Kott, Phillip S. (1999) \emph{"The Extended Delete-A-Group Jackknife"}. Bulletin of the International Statistical Instititute. 52nd Session. Contributed Papers. Book 2, pp. 167-168.

Kott, Phillip S. (2001) \emph{"The Delete-A-Group Jackknife"}. Journal of Official Statistics, Vol.17, No.4, pp. 521-526.
}
\seealso{
\code{\link{desc}} for a concise description of \code{kott.design} objects, \code{\link{kottby}}, \code{\link{kott.ratio}}, \code{\link{kott.regcoef}}, \code{\link{kott.quantile}} and \code{\link{kottby.user}} for calculating estimates and standard errors, \code{\link{kottcalibrate}} for calibrating replicate weights.
}
\examples{
# Creation of kott.design objects starting with survey data sampled
# with different sampling designs (actually the survey data frame is
# always the same: the examples serve the purpose of illustrating
# the syntax).

data(data.examples)

# Two-stage stratified cluster sampling design (notice the presence of
# lonely PSUs):
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~stratum,
      weights=~weight,nrg=15)
desc(kdes)


# The same using collapsed strata (SUPERSTRATUM variable) to remove
# lonely PSUs:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)
desc(kdes)


# Same design, but using the self.rep.str argument to identify
# the SR strata (actually towcod identifies the
# "variance PSUs" by construction):
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15,self.rep.str=~sr)
desc(kdes)


# Two stage cluster sampling (no stratification):
kdes<-kottdesign(data=example,ids=~towcod+famcod,weights=~weight,nrg=15)
desc(kdes)


# One-stage stratified cluster sampling:
kdes<-kottdesign(data=example,ids=~towcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)
desc(kdes)


# Stratified independent sampling design:
kdes<-kottdesign(data=example,ids=~key,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)
desc(kdes)
}
\keyword{survey}