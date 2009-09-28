\name{constrReduce}
\alias{constrReduce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Impose constraints on a QCA solution}
\description{
\code{constrReduce} allows you to include a set of configurations
to the a solution, or exclude a set of configurations from a solution,
then return a new solution.
}
\usage{
constrReduce(object, exclude = NULL, include = NULL,necessary = NULL)
}
\arguments{
  \item{object}{ An object of class "QCA".}
  \item{exclude}{ A data frame, each row represent one configuration.}
  \item{include}{ A data frame, each row represent one configuration.}
  \item{necessary}{ A list, specifying the necessary conditions.}
}
\details{
  Sometime, you may encounter contraditory simplifying assumptions. In
  that case, you may want to exclude the CSAs to attain a more reasonable
  solution.

  In other case, you may attain a solution without including all the
  remainders, latter you may want to include a small number of
  remainders in order to get a intermediate solution.

  Either case you can use constrReduce. Usually, the exclude and include
  argument should repsent a set of configurations. In these case, there
  should be no 'NA' in the data frame. However, it is NOT wrong to
  include NA in the data frame. When there is NA, it means a multiple
  configurations.

  If you attain a solution without including remainders, you can see if
  there is necessary condition. If it does, then you may want to include
  the remainders containning the necessary conditions only. There are
  two ways to do it. First to attain a solutions including
  remainders. then 1) impose constraints by necessary argument. 2) use
  \code{SA} to get the simplyfing solution, find out the SAs not
  containing the necessary conditions and exclude them. See the second
  example.
}
\value{
  An object of class "QCA". It is esentailly a list of 10
  components. See \code{reduce} for more details.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Rongggui HUANG}
%\note{ ~~further notes~~ }
\seealso{ \code{\link{reduce}}}
\examples{
\dontrun{
## example 1
data(Yamasaki,package="QCA")
ans0 <- reduce(Yamasaki,"AGENDA",names(Yamasaki)[1:5],"negative","include") ## 5 solutions
ans1 <- reduce(Yamasaki,"AGENDA",names(Yamasaki)[1:5],"positive","include")
(csa2 <- CSA(ans0[2],ans1)) ## get and show CSAs
(ans02 <- constrReduce(ans0[2],exclude=csa2$solutions[[1]])) ## impose constraint
CSA(ans02,ans1) ## no CSA now

## example 2
data(Osa,package="QCA") ## QCA package is required to run this example
conditions <- c("DYNA","ACCES","INFLU","ELITE","SOCIAL")
a <- reduce(mydata = Osa, outcome = "OUT", conditions = conditions, explain = "negative", remainders = "ex", contradictions = "negative")
b <- reduce(mydata = Osa, outcome = "OUT", conditions = conditions, explain = "negative", remainders = "include", contradictions = "negative")
sa <- SA(b) ## simplifying assumptions
constrReduce(a,inc=sa$solutions[[1]]) ## == b
constrReduce(b,exc=sa$solutions[[1]]) ## == a
constrReduce(b,necess=list(ACCES=0,ELITE=0,social=0)) ## method 1 of imposing necessary conditions
## inspect the SAs (that is sa), 3-6 SAs are against the necessary conditions
constrReduce(b,exclude=sa$solutions[[1]][3:6,]) ## method 2 of imposing necessary conditions
constrReduce(a,include=sa$solutions[[1]][1:2,]) ## the same.
}
}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
