\name{KatzHauMahoney}
\Rdversion{1.1}
\alias{KatzHauMahoney}
\docType{data}
\title{
The "great reversal" data set.
}
\description{
Data set from Kate, Hau and Mahoney (2005) to explain the "great
reversal" in Spanish America.
}
\usage{data(KatzHauMahoney)}
\format{
  A data frame with 14 observations on the following 12 variables.
  \describe{
    \item{\code{Country}}{a factor with levels \code{Argentina}
      \code{Bolivia} \code{Chile} \code{Columbia} \code{CostaRica}
      \code{Ecuador} \code{ElSalvador} \code{Guatemala} \code{Honduras}
      \code{Mexico} \code{Nicaragua} \code{Paraguay} \code{Peru}
      \code{Uruguay} \code{Venezuela}}
    \item{\code{EDF}}{Economically Developed Country. The fuzzy-set
      outcome.}
    \item{\code{EUDF}}{Economically Underdeveloped Country. The
      fuzzy-set outcome}
    \item{\code{GDP}}{GDP in 1990. Statistical outcome}
    \item{\code{SDF}}{Socially Developed Country. The fuzzy-set outcome}
    \item{\code{SUDF}}{Socially Underdeveloped Country. The fuzze-set
      outcome}
    \item{\code{SDI90}}{Socially Development Index in 1990. Statistical
      outcome}
    \item{\code{DIPF}}{Dense Indigenous Population, fuzzy-set score}
    \item{\code{LIEF}}{Labor-Intensive Estates, fuzzy-set score}
    \item{\code{MTEF}}{Mineral/Tropical Exports, fuzzy-set score}
    \item{\code{SLF}}{Strong Liberals, fuzzy-set score}
    \item{\code{SCF}}{Strong Conservatives, fuzzy-set score}
  }
}
\details{
  This study tries to explain the "great reversal",in Spanish
  America. The great reversal refers to the fact that "from 1750
  to 1900, the most marginal colonial territories often became the
  region's wealthiest countries, whereas the most central colonial
  territories often became the region's poorest countries."
  }
\source{
  Katz, A.; vom Hau, M. & Mahoney, J. Explaining the Great Reversal in
  Spanish America: Fuzzy-Set Analysis Versus Regression Analysis
  Sociological Methods and Research, 2005, 33 (4) , 539-573
}
\examples{
data(KatzHauMahoney)
}
\keyword{datasets}
