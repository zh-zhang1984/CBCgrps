\name{dt}
\alias{dt}
\docType{data}
\title{
simulated dataset as a working example
}
\description{
A data frame with 500 observations on the following 8 variables. Missing values are also present
}
\usage{data("dt")}
\format{
  A data frame with 500 observations on the following 8 variables.
  \describe{
    \item{\code{crp}}{a numeric vector, C-reactive protein measured in mg/l, with missing values}
    \item{\code{vaso}}{a factor with two levels \code{Yes} \code{No}, indicating the use of vasopressor or not}
    \item{\code{wbc}}{a numeric vector, white blood cell count}
    \item{\code{lac}}{a numeric vector, serum lactate}
    \item{\code{age}}{a numeric vector, age in years}
    \item{\code{type}}{a factor with levels \code{surgery} \code{medical} \code{emergency}}
    \item{\code{gender}}{a factor with levels \code{female} \code{male}}
    \item{\code{mort}}{an integer with two values \code{1} \code{0}}
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
data(dt)
## maybe str(df) ; plot(df) ...
}
\keyword{datasets}
