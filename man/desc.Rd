\name{desc}
\alias{desc}
\title{Description of replicated objects}
\description{
Concisely describes a \code{kott.design} object.
}
\usage{
desc(deskott, descfun = NULL, ...)
}
\arguments{
  \item{deskott}{Object of class \code{kott.design} containing the replicated survey data.}
  \item{descfun}{Optional description function to be used; must accept a \code{data.frame} object as first argument.}
  \item{\dots}{Additional parameters to be passed to \code{descfun}.}
}
\details{
This function prints a concise description (i) of the sampling design for the original survey data and (ii) of the replication process these data have undergone.

The optional argument \code{descfun} allows to specify an \R function (like \code{head}, \code{str}, \code{summary}, \ldots) to be used to analyse, describe, or summarise the data frame contained in \code{deskott}.
}
\value{
The return value depends on the \code{descfun} parameter. If not specified (the default option), \code{desc} does not return any value.
}
\author{Diego Zardetto}
\examples{
data(data.examples)

# Creation of a kott.design object:
kdes<-kottdesign(data=example,ids=~towcod+famcod,strata=~SUPERSTRATUM,
      weights=~weight,nrg=15)

# Concise description:
desc(kdes)

# Display first rows of kdes data:
desc(kdes,head)

# Ask essential information on kdes internal structure:
desc(kdes,str)


# Creation of a kott.cal.design object:
kdescal04p<-kottcalibrate(deskott=kdes,df.population=pop04p,
            calmodel=~x1+x2+x3-1,partition=~regcod,calfun="logit",
            bounds=bounds,aggregate.stage=2)

# Concise description:
desc(kdescal04p)

# Display first rows of kdescal04p data:
desc(kdescal04p,head)
}
\keyword{survey}