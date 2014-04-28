\name{directCalibration}
\alias{directCalibration}
\title{
  Calibration in fuzzy set QCA.
}
\description{
  Direct method of calibration in fuzzy set QCA.
}
\usage{
directCalibration(x, fullin, fullout, crossover, infz = 0.953, 
                 outfz = 0.047, details = FALSE)
}
\arguments{
  \item{x}{
    an interval variable to be converted to fuzzy set score.
}
  \item{fullin}{
    scalar indicating full membership.
}
  \item{fullout}{
    scalar indicating the full non-membership.
}
  \item{crossover}{
  scalar indicating cross over point.
}
  \item{infz}{
  fuzzy set score of full membership, corresponding to fullin.
}
  \item{outfz}{
  fuzzy set score of full non-membership, corresponding to fullout.
}
  \item{details}{
    logical.
}
}
\details{
  Usually, the sum of infz and outfz should be 1. In Ragin (2008), infz
  is 0.953 and outfz is 0.047, which lead to the log odd of membership
  being (proximitely) 3.0 and -3.0.
}
\value{
  When details is TRUE, it the result is a data frame. Otherwise, the
  value is a numeric vector of fuzzy set score.
}
\references{
  Ragin, Charles C. 2008. "Redesigning social inquiry: fuzzy sets and
  beyond." Chapter 5. Chicago: University of Chicago Press.
}
\author{
  Ronggui HUANG
}
\examples{
 directCalibration(Lipset_fs$Developed,800,300,550)
}
