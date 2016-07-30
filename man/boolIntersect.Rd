\name{boolIntersect}
\alias{boolIntersect}
\title{
Boolean intersection
}
\description{
Boolean intersection
}
\usage{
boolIntersect(..., string = TRUE)
}
\arguments{
  \item{\dots}{
string represented implicants.
}
  \item{string}{
whether returned the string represented outcome or data.frame represented outcome.
}
}
\details{
It is an operator of boolean intersection.
}
\value{
A string or data frame.
}
\author{
Ronggui HUANG}
\examples{
boolIntersect("a*B*D","a*B*c")
}