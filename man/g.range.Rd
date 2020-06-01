\name{g.range}
\alias{g.range}
\title{Range of g-weights}
\description{
Computes the range of the ratios between calibrated weights and direct weights (\emph{g-weights}) for the original sample and all its replicates.
}
\usage{
g.range(cal.deskott)
}
\arguments{
  \item{cal.deskott}{Object of class \code{kott.cal.design}.}
}
\details{
This function computes, for the original sample and all its replicates, the smallest interval which contains the ratios between calibrated weights and direct weights.
}
\value{
A \code{data.frame} object.
}
\author{Diego Zardetto}
\seealso{
\code{\link{kottcalibrate}} for calibrating replicate weights and \code{\link{bounds.hint}} to obtain an hint for calibration problems where range restrictions are imposed on the \emph{g-weights}.
}
\examples{
data(data.examples)

# Creation of a kott.design object:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)

# Calibration:
kdescal04p<-kottcalibrate(deskott=kdes,df.population=pop04p,
            calmodel=~x1+x2+x3-1,partition=~regcod,calfun="logit",
            bounds=bounds,aggregate.stage=2)

# Computation of g-weights range:
g.range(kdescal04p)
}
\keyword{survey}