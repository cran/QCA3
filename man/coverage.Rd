\name{coverage}
\alias{coverage}
\alias{coverage.QCA}
\alias{coverage.default}
\alias{consistency}
\title{
Coverage and consistency.
}
\description{
  Coverage and consistency.
}
\usage{

coverage(x, ...)

\method{coverage}{default}(x, y, alternative = c("less", "greater"),...)

\method{coverage}{QCA}(x,traditional=TRUE,...)

consistency(x, y, alternative = c("less", "greater"))
}
\arguments{
  \item{x}{a numeric vector of fuzzy set score for default method of
    coverage and consistency. A QCA object for QCA method of coverage.
}
\item{y}{
  a numeric vector of fuzzy set score.
}
  \item{alternative}{
    specify the relationship between x and y.
  }
  \item{traditional}{logical, use traditional symbol to present the
    solution when it is TRUE. Otherwise, use Tosmana-style symbol.}
  \item{\ldots}{not used currently.}
}
\details{
  Given two fuzzy set score variables x and y, coverage returns the coverage
  of x by y, given their relationship; and consistency returns  the
  consistency of x by y.
  
  When alternative is "less", it returns coverage(x<=y). This is
  coverage for sufficent test. When "greater", it returns
  coverage(x>=y), it is coverage for necessary condition.

  This value can be used to evaluate the importance or relevance of x as
  a necessary condition for Y. When coverage(x>=y) is small, the
  constraining effect of x on y is negligible. In other words, small
  value of coverage(x>=y) indicates en empirically trivial necessary
  condition (Ragin, 2008: 61).
}
\value{
  For default method of coverage and consistency, the value is a
  length-1 numeric vector.

  For QCA method of coverage, it is a list of data frame. Each data
  frame represent coverage of one solution. It contains three columns,
  namely implicant, raw and unique (coverage in terms of number of case).
}
\references{
  Ragin, C. 2006. Set Relations in Social Research: Evaluating Their
  Consistency and Coverage Political Analysis,14(3):291-310.
  
  Ragin, C. 2008, Redesigning social inquiry: fuzzy set and
  beyond. Chicago: The University of Chicago Press.
}
\author{
  Ronggui HUANG
}
%\seealso{}
\examples{
## evaluate the congerage of necessary conditions for "no democratic consolidation"
## see help(CoD) for more information.
## coverage for necessary test for not-CoD (Schneider, 2009:77)

coverage(1-CoD$econdev,1-CoD$cod,'greater')
# not socioeconomically developed
coverage(1-CoD$demex,1-CoD$cod,'greater')
# no previous democratic experience
coverage(pmin(1-CoD$demex,1-CoD$econdev),1-CoD$cod,'greater')
# not socioeconomic developed society without prior democratic experience

## coverage method for QCA object
tb2<-reduce(policyChange~soc+pol+ind+are+foc,Yamasaki2009,cases="country")
coverage(tb2)
}

