\name{df}
\alias{df}
\docType{data}
\title{
simulated dataset as a working example
}
\description{
A data frame with 1000 observations on the following 7 variables.
}
\usage{data("df")}
\format{
  A data frame with 1000 observations on the following 7 variables.
  \describe{
    \item{\code{crp}}{a numeric vector, C-reactive protein measured in mg/l}
    \item{\code{hb}}{a numeric vector, hemoglobin}
    \item{\code{ddimer}}{a numeric vector}
    \item{\code{wbc}}{a numeric vector, white blood cell}
    \item{\code{comorbid}}{a factor with levels \code{cirrhosis} \code{COPD} \code{diabetes} \code{heartfailure} \code{hypertension} \code{renalfailure} \code{stroke}}
    \item{\code{sex}}{a factor with levels \code{female} \code{male}}
    \item{\code{mort}}{a factor with levels \code{alive} \code{dead}}
  }
}
\details{
The dataset is generated as a working example without clinical relevance. 
}
\source{
Simulated dataset without sources.
}
\references{
simulated dataset without reference.
}
\examples{
data(df)
## maybe str(df) ; plot(df) ...
}
\keyword{datasets}
