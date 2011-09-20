\name{reduce}
\alias{reduce}
\alias{reduce.default}
\alias{reduce.data.frame}
\alias{reduce.truthTable}
\alias{reduce.formula}
\title{ Boolean miniziation for csQCA, mvQCA and fsQCA }
\description{
  This is the core funtion for QCA (Qualitative Comparative
  Analysis). Given the outcome and conditions, it returns an object of
  class 'QCA', which contains all the possible configurations leading to
  the outcome. It can handle various kinds of QCA., namely csQCA,
  mvQCA, fsQCA and csTQCA.
}
\usage{
reduce(x,...)

## default method is an alias of truthTable method.

\method{reduce}{truthTable}(x, explain = c("positive", "negative"),
       remainders = c("exclude","include"),
       contradictions = c("remainders","positive","negative"),
       dontcare = c("remainders", "positive", "negative"),
       keepTruthTable = TRUE,...)

\method{reduce}{data.frame}(x, outcome, conditions,
        explain = c("positive", "negative"), 
        remainders = c("exclude", "include"),
        contradictions = c("remainders", "positive", "negative"),
        dontcare = c("remainders", "positive", "negative"),
        preprocess = c("cs_truthTable", "fs_truthTable",
        "mv_truthTable"), 
       keepTruthTable = TRUE, ...)

\method{reduce}{formula}(x, data, explain = c("positive", "negative"),
      remainders = c("exclude", "include"),
      contradictions = c("remainders", "positive", "negative"),
      dontcare = c("remainders", "positive", "negative"), 
      preprocess = c("cs_truthTable", "fs_truthTable", "mv_truthTable"), 
      keepTruthTable = TRUE, ...)
}
\arguments{
  \item{x}{a R object, it could be a truthTable, data frame or a formula}
  \item{outcome}{a character string to specify the outcome}
  \item{data}{a data frame, which is not optional.}
  \item{conditions}{a character vector to specify the conditions}
  \item{explain}{a character string specifying the cases to be
    explained. Must one of "positive" or "negative".}
  \item{remainders}{a character string specifying how to deal with
    remainders. Must one of "exclude" or "include".}
  \item{contradictions}{a character string specifying how to deal with
    contraditory configurations. Must one of "remainders","positive" or
    "negative"}
  \item{dontcare}{a character string specifying how to deal with
    dontcare cases. Must one of "remainders", "positive", "negative"}
  \item{preprocess}{a character string specifying the function for
    preprocessing data, which turns raw data to a truthTable. Must one
    of \code{cs_truthTable}, \code{fs_truthTable} or \code{mv_truthTable}.}
  \item{keepTruthTable}{logical, when TRUE the returned object keeps
    the truthTable}
  \item{\dots}{ other arguments passed to a function.}
}
\details{
  Outcome is the variable to be explained by the conditions. Conditions
  is explanatory variables that may affect the outcome. It is not
  "independent variable" in statistical sense. Configuration is a
  combination of conditions relevant to a given outcome. Remainders are
  configurations that lack empirical instances. Conraditory
  configuration is a configuration whose outcome value is positive[1]
  for some cases and negative[0] for other cases.

  It is good practices to attain the solutions for both positive outcome with
  and without remainders, and negative outcome with and without
  remainders. If a common necessary condition appears in both solutions
  for positive and negative outcome (without remainders), then such
  necessary condition is a trivial necessary condition(Caramani,
  2009:62). It is not necessary to include trivial necessary condition
  in the final solutions.

  It is good practices to generate and examine a truthTable, then use
  truthTable method of reduce to do the boolean minimization.

  Since version 0.0-3, reduce uses enhanced internal function ereduce1
  (which uses enhanced internal function esubset). It has been tested
  and yields the same result (see tests directory for details).
}
\value{
  An object of class "QCA". It is essentailly a list of 10 components.
  \item{solutions}{ a list of data.frame, each data frame represents one solution.}
  \item{commonSolutions}{A list of lenth nrow(solutionsIDX). For each
  row of solutionsIDX, if the primeImplicants index are the same for all
  solutions, then the index is return. Otherwise, it is NULL.}
  \item{solutionsIDX}{a matrix. Each column represents one
  solution. The number of the matrix is row index of primeImplicants.}
  \item{primeImplicants}{A matrix of prime implicants.}
  \item{truthTable}{a truthTable if keepTruthTable is TRUE, otherwise NULL.}
  \item{explained}{A data frame, representing the configuration of conditions for explained cases. Note it is not on basis of case but basis of configuration.}
  \item{idExclude}{integer vector. id of observed configurations that are
  excluded from minimization. The meaning of id is equivalent to the
  line number of a configuration discussed in Dusa (2007).} 
  \item{nlevels}{a integer vector, the number of levels of each condition.}
  \item{PIChart}{a prime implicants charts, constructed according to
  primeImplicants and explained. It is a logic matrix with dimension of
  nrow(primeImplicants)x ncol(explained). It is TRUE if the
  corresponding primeImplicant covers the corresponding explained.}
  \item{call}{the matched call.}
}
\references{
  Caramani, Daniele. 2009. "Introduction to the comparative method with
  Boolean algebra." Sage.

  Cronqvist, Lasse and  Berg-Schlosser, Dirk. 2009. Multi-Value QCA
  (mvQCA). In Configuraional comparative Methods: qualitative
  comparative analysis (QCA) and related techniques. ed by Benoit RiHoux
  and Charles Ragin. Sage.

  Ragin. Charles. 2009. Qualitative Comparative Analyais Using Fuzzy
  Sets (fsQCA). In Configuraional comparative Methods: qualitative
  comparative analysis (QCA) and related techniques. ed by Benoit RiHoux
  and Charles Ragin. Sage.

  Rihoux, Benoit and De Meur, Gisele. 2009. Crip-Set Qualitative Comparative Analysi
  (csQCA). In Configuraional comparative Methods: qualitative
  comparative analysis (QCA) and related techniques. ed by Benoit RiHoux
  and Charles Ragin. Sage.


  Dusa, Adrian. 2007. Enhancing Quine-McCluskey,
  \url{http://www.compasss.org/files/WPfiles/Dusa2007a.pdf}

  Ragin, Charles. 2000. Fuzzy-Set Social Science. University Of Chicago Press.

  Ragin, Charles. 1987. The Comparative Method. Moving beyond qualitative
  and quantitative strategies. University of California Press.
}
\author{ Ronggui HUANG}
\note{
 With a 2.4 GHz and 4.0GB PC, it takes about about 0.5 minute for 12
  conditions, 2 minutes for 13 conditions, 4.5 minutes for 14
  conditions, and 9 minutes for 15 conditions. It may take a long time
  to get the solution when there are more conditions. You may use
  \code{\link[QCA]{eqmcc}} if speed becomes an issue for
  \code{\link{reduce}}. The disparity is due to the fact that
  \code{eqmcc} eliminates redundant PIs before solving the PIChart
  (Thanks Adrian for pointting it out), but \code{reduce} does not.
  \code{reduce} is a bit greedy in terms of memory usage, for 15 conditions, it
  uses proximately 500 to 600 Mb memory in typical QCA study. I emphasis
  "typical" because the exact scenario also depends on the number of
  observed configurations.
}
\seealso{
\code{\link[QCA]{factorize}}, \code{\link{SA}}, \code{\link{CSA}},
  \code{\link{constrReduce}}
}
\examples{
if (require(QCA)){
data(Osa,package="QCA")
## QCA package is required to run this example
Osa$OUT2 <- Osa$OUT ## OUT is reserved word in QCA3
## the same as examples of QCA:::qmcc
conditions <- c("DYNA","ACCES","INFLU","ELITE","SOCIAL")
reduce(Osa,"OUT2",conditions,explain="positive",remaind="exclude")
reduce(Osa,"OUT2",conditions,explain="positive",contradictions="positive",remaind="include")
ans <-
  reduce(Osa,"OUT2",conditions,explain="positive",contradictions="negative",remaind="include")
simplifyingAssumption(ans) ## or SA(ans)
reduce(Osa,"OUT2",conditions,explain="negative",contradictions="negative",remaind="include")

## Results of Osa and Corduneanu-Huci (2003)
reduce(Osa,"OUT2",conditions,explain="pos",contrad="neg",remaind="exclude")
# table 1 in page 617
reduce(Osa,"OUT2",conditions,explain="neg",contrad="pos",remaind="exclude")
# table 2 of page 621
reduce(Osa,"OUT2",conditions,explain="positive",contradictions="pos",remaind="incl")
# maximum reduction in page 623
reduce(Osa,"OUT2",conditions[1:4],explain="pos",contradictions="neg",remaind="excl")
# Appendix 2 in page 629
}

## csQCA, mvQCA and fsQCA examples from "Configuraional comparative Methods"
## csQCA
conditions <- c("GNPCAP", "URBANIZA", "LITERACY", "INDLAB", "GOVSTAB")
reduce(Lipset_cs,"SURVIVAL",conditions,explain="positive",remainder="exclude",case="CASEID")
## or use formula
reduce(SURVIVAL~GNPCAP+URBANIZA+LITERACY+INDLAB+GOVSTAB,Lipset_cs,
       explain="positive",remainder="exclude",case="CASEID")
## Formula 1 in Rihoux and De Meur(2009:57)
reduce(Lipset_cs,"SURVIVAL",conditions,explain="negative",remainder="exclude",case="CASEID")
## Formula 3 in Rihoux and De Meur(2009:59)
ans1 <-
 reduce(Lipset_cs,"SURVIVAL",conditions,explain="positive",remainder="include",case="CASEID")
print(ans1) ## Formula 4 in Rihoux and De Meur(2009:60)
SA(ans1) ## 5 simplifying assumptions in p61
ans0 <-
 reduce(Lipset_cs,"SURVIVAL",conditions,explain="negative",remainder="include",case="CASEID")
print(ans0) ## Formula 5 in Rihoux and De Meur(2009:61)
SA(ans0) ## 18 simplifying assumptions

## mvQCA
conditions <- c("GNPCAP", "URBANIZA", "LITERACY", "INDLAB")
reduce(Lipset_mv,"SURVIVAL",conditions,explain="positive",remainder="exclude",case="CASEID",prep="mv_truthTable")
## formula 1 Cronqvist and Berg-Schlosser(2009:80)
ans1 <-
  reduce(Lipset_mv,"SURVIVAL",conditions,explain="positive",remainder="include",case="CASEID",prep="mv_truthTable")
print(ans1) ## formula 2 in Cronqvist and Berg-Schlosser(2009:81)
SA(ans1) ## 9 SAs (see end note 7)
reduce(Lipset_mv,"SURVIVAL",conditions,explain="negative",remainder="exclude",case="CASEID",prep="mv_truthTable")
## formula 3 in Cronqvist and Berg-Schlosser(2009:81)
ans0 <-
  reduce(Lipset_mv,"SURVIVAL",conditions,explain="negative",remainder="include",contrad="positive",case="CASEID",prep="mv_truthTable")
print(ans0) ## formula 4 in Cronqvist and Berg-Schlosser(2009:81)
SA(ans0) ## 7 SAs (see end note 9)

## fsQCA
conditions <- c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ")
reduce(Lipset_fs,"Survived.FZ",conditions,explain="positive",remaind="exclude",prepro="fs",consistency=0.7)
## Formula 1 in Ragin (2009:112)
reduce(Lipset_fs,"Survived.FZ",conditions,explain="positive",remaind="include",prepro="fs",consistency=0.7)
## Formula 2 in Ragin (2009:114)
reduce(Lipset_fs,"Survived.FZ",conditions,explain="negative",remaind="exclude",prepro="fs",consistency=0.7)
## Formula 5 in Ragin (2009:115)
reduce(Lipset_fs,"Survived.FZ",conditions,explain="negative",remaind="include",prepro="fs",consistency=0.7)
## Formula 6 in Ragin (2009:117)
}
