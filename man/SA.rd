\name{simplifyingAssumption}
\alias{simplifyingAssumption}
\alias{SA}
\alias{CSA}
\title{Simplyfing assumptions and contraditory simplifying assumptions}
\description{
  Return the simplifying assumptions and contradictory simplifying assumptions.
}
\usage{
simplifyingAssumption(object, ...) 
SA(object, ...) ## shortcut of simplifyingAssumption
CSA(object1, object0)
}

\arguments{
  \item{object}{An object of class "QCA", which is return from \code{reduce}}
  \item{\dots}{Not used currently}
  \item{object1}{An object of class "QCA" with one solution.}
  \item{object0}{An object of class "QCA" with one solution.}
}
\details{
  Simplyfying assumption is assumption made on the outcome value of a
  logical remainder, so it can be included in the minimization
  procedure. Thus, it is meaning to use \code{SA} and \code{CSA} when
  the object is return by a call to \code{reduce} with remainder
  arugment set to "include".
  
  A contraditory simplifying assumption (CSA) occurs when the same logical
  remainder is used both in the minization of the positive outcome
  configurations and in the minization of the negative outcome
  configuration. The CSA should be solved. An overly heavy presence of
  CSAs is one indicator of problem in the selection of conditions.

  If can object of class "QCA" have multiple solutions, you can use
  \code{[} to extract one of solution, then pass it to \code{CSA}. see
  example section for an example. For object1 and object0, one is the
  solution for explaination of positive case and the other is the
  solution for explaination of negative case.
}
\value{
For \code{SA} and \code{CSA}, the value is an object of class
  c("SA","QCA"). It is a list of 8 components.
}
\references{
  Yamasaki and Rihoux. 2009. A commented review of applications. In
  Configuraional comparative Methods: qualitative comparative analysis
  (QCA) and related techniques. ed by Benoit RiHoux and Charles
  Ragin. Sage.
}
\author{Ronggui HUANG}
\seealso{\code{\link{reduce}} \code{\link{constrReduce}}}
\examples{
\dontrun{
data(Yamasaki,package="QCA")
cond <- names(Yamasaki)[1:5]
ans0 <- reduce(Yamasaki,"AGENDA",cond,"negative","include") ## 5 solutions 
ans1 <- reduce(Yamasaki,"AGENDA",cond,"positive","include") ## 1 solutions
SA(ans0)
SA(ans1)
CSA(ans0[1],ans1) ## no CSA, please note the subset operation
CSA(ans0[2],ans1) ## ans0[2]-ans0[5] have CSA
}
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
