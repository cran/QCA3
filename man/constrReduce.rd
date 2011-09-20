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

 \code{excludeCSA} conducts Boolean minimization including other
 remainders except those in condictory simplifying assumptions.

 The difference between \code{constrReduce} and \code{excludeCSA} mainly
 lies in how to deal with other remainders when imposing constraints on a QCA
 solution. \code{constrReduce} does not include further remainders,
 while \code{excludeCSA} does.
 
  Sometime, you may encounter contraditory simplifying assumptions. In
  that case, you may want to exclude the CSAs to attain a more reasonable
  solution. In this case, \code{excludeCSA} is the most suitable way to
  go, which can make QCA easier. However, it does not guarantee a final
  solution, in particular when there are multiple solutions to both
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
  esentailly a list of10 components. See \code{\link{reduce}} for more
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

if (require(QCA)){
##Use easy counterfactuals to get an intermediate solution (Ragin 2008:chapter 9)
data(Stokke,package="QCA")
comp <- reduce(Stokke,"Y",c("A","C","S","I","R"),explain="positive")
pars <- reduce(Stokke,"Y",c("A","C","S","I","R"),explain="positive",remaind="include")
sa <- SA(pars)
## determins easy counterfactuals
idx <- c(12,14,## A*S*R must be kept
         14, ## A*C*S*i must be kept
         3,5,8,10,12,14 ## A*i must be kept
         )
idx <- unique(idx) ## index of easy counterfactuals
easy <- sa$solution[[1]] [idx,]
difficult <- sa$solution[[1]] [-idx,]
constrReduce(comp,include=easy)
constrReduce(pars,exclude=difficult) ## the last two are equivalent.

## another way is to manually construct the easy counterfactuals
easy2 <- rbind(
c(1,-9,1,-9,1), # A*S*R
c(1,1,1,0,-9), # A*C*S*i
c(1,-9,-9,0,-9) # A*i
)
easy2 <- as.data.frame(easy2)
constrReduce(comp,include=easy2)
## end of Ragin (2009:chapter 9) example

## example 1
data(Yamasaki,package="QCA")
ans0 <- reduce(Yamasaki,"AGENDA",names(Yamasaki)[1:5],"negative","include") ## 5 solutions
ans1 <- reduce(Yamasaki,"AGENDA",names(Yamasaki)[1:5],"positive","include")
(csa2 <- CSA(ans0[2],ans1)) ## get and show CSAs
(ans02 <- constrReduce(ans0[2],exclude=csa2$solutions[[1]]))
## impose constraint
CSA(ans02,ans1) ## no CSA now

ans02b <- excludeCSA(ans0[2],csa2)
CSA(ans02b,ans1) ## no CSA
## note that ans02b is more parsimonious than ans02
## since ans02b includes further reminders.

## example 2
data(Osa,package="QCA") ## QCA package is required to run this example
conditions <- c("DYNA","ACCES","INFLU","ELITE","SOCIAL")
Osa$OUT2 <- Osa$OUT ## OUT is reserved word in QCA3.
a <- reduce(Osa, outcome = "OUT2", conditions = conditions, explain = "negative", remainders = "ex", contradictions = "negative")
b <- reduce(Osa, outcome = "OUT2", conditions = conditions, explain = "negative", remainders = "include", contradictions = "negative")

# there are two solutions in b, let's focus on the first only.
b1 <- b[1]
sa <- SA(b1) ## simplifying assumptions
constrReduce(a,inc=sa$solutions[[1]]) ## == b1
constrReduce(b1,exc=sa$solutions[[1]]) ## ==a1
}
}

