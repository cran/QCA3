\name{reduce}
\alias{reduce}
\alias{reduce.truthTable}
\alias{reduce.default}
\alias{print.QCA}
\alias{summary.QCA}
\alias{[.QCA}
\title{ Boolean miniziation for csQCA, mvQCA and fsQCA }
\description{
  This is the core funtion for QCA (Qualitative Comparative
  Analysis). Given the outcome and conditions, it returns an object of
  class PI, which contains all the possible configurations leading to
  the outcome. It can handle all three kinds of QCA., namely csQCA,
  mvQCA, and fsQCA.
}
\usage{
reduce(mydata,...)

\method{reduce}{truthTable}(mydata, explain = c("positive", "negative"),
       remainders = c("exclude","include"),
       contradictions = c("remainders","positive","negative"),
       dontcare = c("remainders", "positive", "negative"), 
       keepTruthTable = TRUE,...)

\method{reduce}{default}(mydata, outcome, conditions,
       explain = c("positive", "negative"),
       remainders = c("exclude", "include"),
       contradictions = c("remainders","positive", "negative"),
       dontcare = c("remainders", "positive", "negative"),
       preprocess = c("cs_truthTable","fs_truthTable", "pass"),
       nlevels = rep(2, length(conditions)), keepTruthTable = TRUE, ...)

\method{print}{QCA}(x, traditional = TRUE, show.truthTable = TRUE, ...)

\method{summary}{QCA}(object, traditional = TRUE, show.case = TRUE, ...)

\method{[}{QCA}(object, which)
}
\arguments{
  \item{mydata}{a data frame}
  \item{outcome}{a character string to specify the outcome}
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
    of \code{cs_truthTable}, \code{fs_truthTable} or \code{pass}}
  \item{nlevels}{a integer vector, specifying number of levels for the
    corresponding conditions. For csQCA and fsQCA, it is always
    \code{rep(2,length(conditions))}}
  \item{keepTruthTable}{logical, when TRUE the returned object keeps
    the truthTable}
  \item{\dots}{ For \code{reduce}, arguements passed to preprocess
    function; for \code{print.QCA} and \code{summary.QCA}, currently not
    used.}
  \item{x}{an object of class 'QCA', which is usually returned from
    \code{reduce}.}
  \item{traditional}{logical, use traditional symbol when it is
    TRUE. Otherwise, use Tosmana-style symbol.}
  \item{show.truthTable}{logical, show truthTable when it is TRUE. Of
    course, it has effect only when the 'keepTruthTable' argument of
    \code{reduce} is set to TRUE.}
  \item{object}{an object of class 'QCA', which is usually returned from
    \code{reduce}.}
  \item{show.case}{logical, show case names when it is TRUE.}
  \item{which}{numeric vector, indices specifying elements to
extract. Extraction of a solution or (prime implicant) is essentially a
extraction on a list. you can refer to \code{[} for more details.}
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

  The traditional way uses upper-case letters representing 1 and and
  lower-case letters reprensenting 0. The Tosmana-style uses
  \code{condition{value}} to represent the prime implicants.
}
\value{
  An object of class "QCA". It is essentailly a list of 10 components.
}
\references{
  Caramani, Daniele 2009. "Introduction to the comparative method with
  Boolean algebra." SAGE.
  
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

  
  Dusa, Adrian 2007 Enhancing Quine-McCluskey,
  \url{http://www.compasss.org/files/wpfiles/Dusa2007a.pdf}
  
  Ragin, Charles. 2000. Fuzzy-Set Social Science. University Of Chicago Press.

  Ragin, Charles. 1987. The Comparative Method. Moving beyond qualitative
  and quantitative strategies. University of California Press.
}
\author{ Ronggui HUANG}
\note{
It takes about 6 minutes for 12 conditions. It may take a long time to
get the solution when there are more conditions. You may use
  \code{\link[QCA]{eqmcc}} if speed becomes an issue for \code{\link{reduce}}.
}
\seealso{  \code{\link[QCA]{factorize}}, \code{\link{SA}},
  \code{\link{CSA}}, \code{\link{constrReduce}}}
\examples{
\dontrun{
data(Osa,package="QCA") ## QCA package is required to run this example
## the same as examples of QCA:::qmcc
conditions <- c("DYNA","ACCES","INFLU","ELITE","SOCIAL")
reduce(Osa,"OUT",conditions,explain="positive",remaind="exclude")
reduce(Osa,"OUT",conditions,explain="positive",contradictions="positive",remaind="include")
ans <- reduce(Osa,"OUT",conditions,explain="positive",contradictions="negative",remaind="include")
simplifyingAssumption(ans) ## or SA(ans)
reduce(Osa,"OUT",conditions,explain="negative",contradictions="negative",remaind="include")

## Results of Osa and Corduneanu-Huci (2003)
reduce(Osa,"OUT",conditions,explain="pos",contrad="neg",remaind="exclude") # table 1 in page 617
reduce(Osa,"OUT",conditions,explain="neg",contrad="pos",remaind="exclude") # table 2 of page 621
reduce(Osa,"OUT",conditions,explain="positive",contradictions="pos",remaind="incl") # maximum reduction in page 623
reduce(Osa,"OUT",conditions[1:4],explain="pos",contradictions="neg",remaind="excl") # Appendix 2 in page 629
}

## csQCA, mvQCA and fsQCA examples are from "Configuraional comparative Methods"
## csQCA
conditions <- c("GNPCAP", "URBANIZA", "LITERACY", "INDLAB", "GOVSTAB")
reduce(Lipset_cs,"SURVIVAL",conditions,explain="positive",remainder="exclude",case="CASEID")
## Formula 1 in Rihoux and De Meur(2009:57)
reduce(Lipset_cs,"SURVIVAL",conditions,explain="negative",remainder="exclude",case="CASEID")
## Formula 3 in Rihoux and De Meur(2009:59)
ans1 <- reduce(Lipset_cs,"SURVIVAL",conditions,explain="positive",remainder="include",case="CASEID")
print(ans1)
## Formula 4 in Rihoux and De Meur(2009:60)
SA(ans1) ## 5 simplifying assumptions in p61
ans0 <- reduce(Lipset_cs,"SURVIVAL",conditions,explain="negative",remainder="include",case="CASEID")
print(ans0)
## Formula 5 in Rihoux and De Meur(2009:61)
SA(ans0) ## 18 simplifying assumptions

## mvQCA
conditions <- c("GNPCAP", "URBANIZA", "LITERACY", "INDLAB")
reduce(Lipset_mv,"SURVIVAL",conditions,explain="positive",remainder="exclude",case="CASEID",nlevels=c(3,2,2,2))
## formula 1 Cronqvist and Berg-Schlosser(2009:80)
ans1 <-
reduce(Lipset_mv,"SURVIVAL",conditions,explain="positive",remainder="include",case="CASEID",nlevels=c(3,2,2,2))
print(ans1)
## formula 2 in Cronqvist and Berg-Schlosser(2009:81)
SA(ans1) ## 9 SAs (see end note 7)
reduce(Lipset_mv,"SURVIVAL",conditions,explain="negative",remainder="exclude",case="CASEID",nlevels=c(3,2,2,2))
## formula 3 in Cronqvist and Berg-Schlosser(2009:81)
ans0 <-
reduce(Lipset_mv,"SURVIVAL",conditions,explain="negative",remainder="include",contrad="positive",case="CASEID",nlevels=c(3,2,2,2))
print(ans0)
## formula 4 in Cronqvist and Berg-Schlosser(2009:81)
SA(ans0) ## 7 SAs (see end note 9)

## fsQCA
conditions <- c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ")
reduce(mydata=Lipset_fs,"Survived.FZ",conditions,explain="positive",remaind="exclude",prepro="fs",consistency=0.7)
## Formula 1 in Ragin (2009:112)
reduce(mydata=Lipset_fs,"Survived.FZ",conditions,explain="positive",remaind="include",prepro="fs",consistency=0.7)
## Formula 2 in Ragin (2009:114)
reduce(mydata=Lipset_fs,"Survived.FZ",conditions,explain="negative",remaind="exclude",prepro="fs",consistency=0.7)
## Formula 5 in Ragin (2009:115)
reduce(mydata=Lipset_fs,"Survived.FZ",conditions,explain="negative",remaind="include",prepro="fs",consistency=0.7)
## Formula 6 in Ragin (2009:117)
}
%\keyword{ ~kwd1 }
