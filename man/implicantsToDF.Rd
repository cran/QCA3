\name{implicantsToDF}
\alias{implicantsToDF}
\title{
Conversion of implicants to a Data Frame representation
}
\description{
Conversion of implicants to a Data Frame representation
}
\usage{
implicantsToDF(x, conditions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A string presenting implicants.
}
  \item{conditions}{
  A character vector of all condition names.
}
}
\value{
A data frame, whose columns are defined by conditions.
}
\author{
Ronggui HUANG
}
\seealso{
\link{constrReduce}
}
\examples{
implicantsToDF(x="A*S*R+A*C*S*i+A*i",conditions=c("A","C","S","I","R"))
}
