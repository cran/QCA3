\name{Lipsets}
%\alias{Lipsets}
\alias{Lipset}
\alias{Lipset_fs}
\alias{Lipset_cs}
\alias{Lipset_mv}
\docType{data}
\title{Breakdown/survival of democracy in inter-war Europe}
\description{
This is the fuzzy set version of Lipset data on Breakdown/survival of democracy in inter-war Europe.
}
%\usage{data(Lipset_fs)}
\format{
  A data frame with 18 observations on the following 13 variables.
  \describe{
    \item{\code{Country}}{Name of country}
    \item{\code{Survived}}{If a country survives the economic and political upheavals of this period.}
    \item{\code{Survived.FZ}}{Fuzzy set score of Survived}
    \item{\code{Developed}}{The degree of development.}
    \item{\code{Developed.FZ}}{Fuzzy set score of Developed.}
    \item{\code{Urban}}{The degree of urbanization.}
    \item{\code{Urban.FZ}}{Fuzzy set score of Urban.}
    \item{\code{Literate}}{The degree of literate a county is.}
    \item{\code{Literate.FZ}}{Fuzze set score of Literate.}
    \item{\code{Industrial}}{The degree of industrialization.}
    \item{\code{Industrial.FZ}}{Fuzze set score of Industrial}
    \item{\code{Unstable}}{The degree of political instablity.}
    \item{\code{Stable.FZ}}{Fuzzy set score of political stability.}
  }
}
\details{
 The data set is from Ragin(2009:95), the details about the dataset is
 described in page 93.
}
\source{
  Ragin. Charles. 2009. Qualitative Comparative Analyais Using
Fuzzy Sets (fsQCA). In Configuraional comparative Methods: qualitative
comparative analysis (QCA) and related techniques. ed by Benoit RiHoux
and Charles Ragin. Sage.
}
\examples{
conditions <- c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ")
reduce(Lipset_fs,"Survived.FZ",conditions,explain="positive",remaind="exclude",prepro="fs",consistency=0.7)
## Formula 1 in page 112
reduce(Lipset_fs,"Survived.FZ",conditions,explain="positive",remaind="include",prepro="fs",consistency=0.7)
## Formula 2 in page 114
reduce(Lipset_fs,"Survived.FZ",conditions,explain="negative",remaind="exclude",prepro="fs",consistency=0.7)
## Formula 5 in page 115
reduce(Lipset_fs,"Survived.FZ",conditions,explain="negative",remaind="include",prepro="fs",consistency=0.7)
## Formula 6 in page 117
}
\keyword{datasets}
