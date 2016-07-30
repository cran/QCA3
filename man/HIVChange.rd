\name{HIVChange}
\alias{HIVChange}
\docType{data}
\title{data about HIV prevalence in Sub-Saharan Africa}
\description{
Data set from Cronqvist and Berg-Schlosser(2006), examining the HIV
prevalence in Sub-Saharan Africa.
}
\usage{data(HIVChange)}
\format{
  A data frame with 12 observations on the following 7 variables.
  \describe{
    \item{\code{LIT00}}{Dichotomized literacy rate of 2000. 1 if above
      50\%.}
    \item{\code{GENDEREQ}}{Dichotomized gender equality index. 1 is
      above 40 index point.}
    \item{\code{MORTALITY}}{Trichotomized cummulated HIV morality rate up
    to 1997, with thresholds of 2\% and 4\%.}
    \item{\code{AGRARGDP}}{Dichotomized variable of share of agarian
      production of GDP. 1 is above 25\%.}
    \item{\code{HIVChange}}{The outcome. "1" if the prevalence rate
      increases. Otherwise it is "0". "C" is indicates contraditory
      configurations.} 
    \item{\code{Country}}{Names of cases}
    \item{\code{NCase}}{Number of cases}
  }
}
\details{
  \code{LIT00} measures socio-economic factors. \code{GENDEREQ} measures
  the situation of women. \code{MORTALITY} measures the awareness of HIV
  threat. \code{AGRARGDP} measures the impact of migration.
}
\references{
  Cronqvist. L. and Berg-Schlosser, D. 2006. Determining The Conditions
  Of Hiv/Aids Prevalence In Sub-Saharan Africa: Employing New Tools Of
  Macro-Qualitative Analysis. In Innovative Comparative Methods For
  Policy Analysis: Beyond The Quantitative-Qualitative Divide. Benoit
  Rihoux and Heike Grimm (Eds).Springer.
}
\examples{
 ## manually construct a truthTable for mvQCA
 conditions <- c("LIT00", "GENDEREQ", "MORTALITY", "AGRARGDP")
 HIVChange$OUT <- HIVChange$HIVChange
 HTT <- HIVChange[,c(conditions,"OUT","NCase","Country")]
 names(HTT) <- c(conditions,"OUT","NCase","Cases")
 ## optional
 mvTT <- list(truthTable=HTT,nlevels=c(2,2,3,2),conditions=conditions)
 ## only some components are rquired
 class(mvTT) <- "truthTable"

 ## Example in p161: not exactly the same solution.
 ##    This one is correct too (tell me if you don't think so)
 reduce(mvTT,expl="pos",remainder="include",contr="negative")
 reduce(mvTT,expl="neg",remainder="include",contr="positive")

 ## Example in p163
 reduce(mvTT,expl="pos",remainder="include",contr="negative")
 reduce(mvTT,expl="neg",remainder="include",contr="negative")
 ## C.A.R is positive, all other three are nagetive
}
\keyword{datasets}
