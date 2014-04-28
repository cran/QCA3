\name{thresholdssetter}
\alias{thresholdssetter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recode continuous variable into binary/multi-value variable}
\description{
If the conditions is continuous variable, it has to be recoded into
binary variable in order to conduct csQCA, or multi-value variable in
order to conduct mvQCA. \code{thresholdsster} helps to recode the
variable using cluster analysis (the default method is 'average linkage').
}
\usage{
thresholdssetter(x, nthreshold = 1,method ="average",dismethod="euclidean",
                 thresholds=NULL, plot=TRUE)
}
\arguments{
  \item{x}{A continouse variable.}
  \item{nthreshold}{The number of thresholds.}
  %\item{value}{logical, returns the recoded variable when TRUE, otherwise
  %returns the thresholods only.}
  \item{method}{method of cluster analysis. See \code{hclust} for more
    details.}
  \item{dismethod}{method argument for \code{dist}.}
  \item{thresholds}{numeric vector of the threshold values. If it is not
  NULL, then recode the variable according to thresholds rather than
  cluster analysis. The threshold values are included in the smaller group.}
  \item{plot}{logical. If TRUE, a bar plot is produced to aid the
    examination of distribution of the recoded variable.}
}
\details{
  In order to dichotomize the data, 1 threshold is needed; to
  trichotomize the data, 2 thresholds are needed; and so on. Yet, it
  should not be a large number, otherwise, limited diversity is a
  problem. It is suggested to be 1 to 3.

  Whenever possible, use theoretically meaningful thresholds. Use
  statistical means only when necessary. The thresholds should make
  theoretical sense, thus you should examine the threshold values before
  going on. Good thresholod values should NOT create very differently
  sized subgroups. Mean and median should be avoided as well.
}
\value{
A numeric sequence from 0 to (nthreshold+1), with the attribute of thresholds.
}
\references{
  Cronqvist, Lasse and  Berg-Schlosser, Dirk. 2009. Multi-Value QCA
  (mvQCA). In Configuraional comparative Methods: qualitative
  comparative analysis (QCA) and related techniques. ed by Benoit RiHoux
  and Charles Ragin. Sage.
    
  Cronqvist, L. 2007. Tosmana user manual. \url{http://www.tosmana.net/tosmana_manual1_3beta.pdf}
}
\author{Ronggui HUANG}
\examples{
thresholdssetter(Lipset$GNPCAP,2)
}

