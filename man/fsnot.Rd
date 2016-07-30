\name{fsnot}
\alias{fsnot}
\alias{fsand}
\alias{fsor}
\title{
  The not opertor for a fuzzy set score.
}
\description{
  The not opertor for a fuzzy set score.
}
\usage{
fsnot(x)

fsand(...)

fsor(...)
}
\arguments{
  \item{x}{
    a variable presenting fuzzy set score
  }
  \item{\dots}{fuzyy set score variables}
}
\details{
  The not operator is define as 1 - x.
}
\value{
a numeric vector of the new fuzzy set score.
}
\author{
Ronggui HUANG
}
\examples{
x = c(0.1,0.5,0.6,0.8)
fsnot(x)

## can be used with fsplot
fsplot(cod~fsnot(econdev),data=CoD)
}