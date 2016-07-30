\name{constrReduce}
\alias{constrReduce}
\alias{excludeCSA}
\title{ Impose constraints on a QCA solution}
\description{
To impose constraints on a QCA solution and returns a new QCA solution.
}
\usage{
constrReduce(object, exclude = NULL, include = NULL,necessary = NULL)

excludeCSA(object, csa) 
}
\arguments{
  \item{object}{ An object of class "QCA".}
  \item{exclude}{ A data frame, each row represent one configuration.}
  \item{include}{ A data frame, each row represent one configuration.}
  \item{necessary}{ A list, specifying the necessary conditions.}
  \item{csa}{an object returned by \code{CSA}.}
}
\details{ 
  \code{constrReduce} allows you to include a set of configurations
  to the a solution, or exclude a set of configurations from a solution,
  then return a new solution.

 \code{excludeCSA} conducts Boolean minimization by freely including other
 remainders except those in condictory simplifying assumptions and
 thoese excluded in the first place (idExclude component of the object).

 The difference between \code{constrReduce} and \code{excludeCSA} mainly
 lies in how to deal with other remainders when imposing constraints on a QCA
 solution. \code{constrReduce} does not include further remainders,
 while \code{excludeCSA} does.
 
  Sometime, you may encounter contraditory simplifying assumptions. In
  that case, you may want to exclude the CSAs to attain a more reasonable
  solution. In this case, \code{excludeCSA} is the most suitable way to
  go, which can make QCA easier. However, it does not guarantee a single
  final solution, in particular when there are multiple solutions to both
  positive and negative outcome. 

  Sometimes, you may attain a solution without including all the
  remainders, latter you may want to include a small number of
  remainders in order to get a intermediate solution. You may include
  raminders in Boolean minimization, but some of the simplifying
  assumptions (see \code{\link{SA}}) are not feasible. You might want to
  exclude these unjustified simplifying assumptions without search for
  further simplifying assumptions. If you attain a solution without
  including remainders, you may want to include easy counterfactuals
  (justifiable remainders) to attain an intermediate solution. In either
  case, \code{constrReduce} is the most suitable way to go. See the
  example section on Ragin (2008: chapter 9) on this usage.

  Usually, the exclude and include argument should repsent a set of
  configurations. In these case, there should be no '-9' in the data
  frame. However, it is NOT wrong to include -9 in the data frame. When
  there is -9, it means a multiple configurations because -9 denotes
  "dont' care".

  % If you attain a solution without including remainders, you can see if
  % there is necessary condition. If it does, then you may want to include
  % the remainders containning the necessary conditions only. The has
  % something to do with the necessary argument in constrReduce. However, it
  % is not a fully-fledged function yet. It needs improvement in future
  % version.
}
\value{
  For \code{constrReduce}, it is an object of class "QCA". It is
  esentailly a list of 11 components. See \code{\link{reduce}} for more
  details. The only difference is the call componenet.
  
}
\author{ Rongggui HUANG}
\references{
  Ragin, Charles C. 2008. "Redesigning social inquiry: fuzzy sets and
  beyond." Chapter 9. Chicago: University of Chicago Press.
}
\seealso{
  \code{\link{reduce}} and  \code{\link{CSA}}
}
\examples{
example(HuangGui2009)
newSol <- excludeCSA(ans2[2],CSA(ans1,ans2[2]))
## ans2 has 3 solutions, only the 3rd has not CSA.
identical(newSol$solutions, ans2[3]$solutions)
## they are the same.

## Use easy counterfactuals to get an intermediate solution (Ragin 2008:chapter 9)
comp <- reduce(Stokke,"success",c("advice","commitment","shadow","inconvenience",
              "reverberation"),explain="positive")
pars <- reduce(Stokke,"success",c("advice","commitment","shadow","inconvenience",
              "reverberation"),explain="positive",remaind="include")
sa <- SA(pars)
## determins easy counterfactuals
## method 1 is to manually construct the easy counterfactuals
easy1 <- rbind(
c(1,-9,1,-9,1), # ADVICE*SHADOW*REVERBERATION
c(1,1,1,0,-9), # ADVICE*COMMITMENT*SHADOW*inconvenience
c(1,-9,-9,0,-9) # ADVICE*inconvenience
)
easy1 <- as.data.frame(easy1)
constrReduce(comp,include=easy1)
## method 2 uses implicantsToDF faciliates such construction
\dontrun{
imp <- "ADVICE*SHADOW*REVERBERATION+
ADVICE*COMMITMENT*SHADOW*inconvenience+
ADVICE*inconvenience"
easy2 <- implicantsToDF(imp, 
          conditions=c("advice","commitment","shadow",'inconvenience',"reverberation"))
constrReduce(comp,include=easy2)
}
## end of Ragin (2009:chapter 9) example
}
