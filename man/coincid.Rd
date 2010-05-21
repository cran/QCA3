\name{coincid}
\alias{coincid}
\alias{print.coincid}
\title{Produces a coincidence matrix}
\description{
Produce a coincidence matrix from a crip set or fuzzy set data.
}
\usage{
coincid(x, standardized = FALSE, use = c("complete", "pairwise"))
%# print method for class coincid
\method{print}{coincid}(x, digits = 3, ...)
}

\arguments{
  \item{x}{For coincide, x is a dataframe of crip set or fuzzy set data, which ranges from
  0 to 1.}
  \item{standardized}{When TRUE, produce the standardized coincidence scores.}
  \item{use}{The method of handling missing data. "complete" means
    listwise deletion and "pairwise" means pairwise deletion.}
  \item{digits}{a non-null value for 'digits' specifies the minimum
    number of significant digits to be printed in values.}
  \item{\ldots}{arguments passed to default method of print.}
}

\details{Coincidence measures the amount of overlap or coincidence between two sets or configurations. The formula is:
  \code{sum(pmin(x,y))/sum(pmax(x,y))}. When standardized is TRUE, the formular is: \code{sum(pmin(x,y))/min(sum(x),sum(y))}.
}
\value{
  A coincidence matrix.
}
\references{
  Ragin, Charles C. 2006. "Set Relations in Social Research: Evaluating Their Consistency and Coverage." Political Analysis 14 (3) : 291-310
}

\author{Ronggui HUANG}
\seealso{\code{\link[QCA3]{suffnec}}}

%\examples{
%}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
