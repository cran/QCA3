\name{fs_tt}
\alias{fs_tt}
\title{Construction of truthTable from fuzzy set score}
\description{
Constructing a truthTable from fuzzy set score. This truthTable is an
object of class tt, which can be passed to procedure \code{\link[QCA]{eqmcc}} to get
the boolean minimization.
}
\usage{
fs_tt(mydata, outcome = "", conditions = c(""),ncases_cutoff = 1,
      consistency_cutoff = 0.8, complete= FALSE, show.cases = TRUE,
      quiet =FALSE, cases =NULL)
}

\arguments{
  \item{mydata}{ A fuzzy set score dataset. All the scores must range
  from 0 to 1.}
  \item{outcome}{ the name of the outcome variable in the dataset.}
  \item{conditions}{the name of the conditions from the dataset.}
% \item{membership_cutoff}{Cutoff point determining the number of case in each
%  corner of the fuzzy set space. Only cases with memebership score
%  greater than membership_cutoff belongs to the combination of conditions (the
%  corner of fuzzy set space).}
  \item{ncases_cutoff}{When number of case is less than cutoff, it will
  be regarded as dontcare configuration.}
  \item{consistency_cutoff}{Cutoff point of consistenty score, cases with
  consistency score greater than cutoff point are regarded as OUT=1. }
  \item{complete}{prints the complete truth table, including the missing combinations}
  \item{show.cases}{show the rownames from the original dataset for each combination of conditions.}
  \item{quiet}{Not used currently.}
  \item{cases}{ name of cases from dataset. If it is NULL and show.cases
  is TRUE, name of cases are derived from row names of dataset.}
}
\details{
  See \code{\link{fs_truthTable}} for more details.
}
\value{
  An object of class "truthTable". refer to \code{\link[QCA]{truthTable}} for more details.
}
\references{
  Ragin. Charles. 2009. Qualitative Comparative Analyais Using
Fuzzy Sets (fsQCA). In Configuraional comparative Methods: qualitative
comparative analysis (QCA) and related techniques. ed by Benoit RiHoux
and Charles Ragin. Sage.
}
\author{ Ronggui HUANG}
\seealso{\code{\link{fs_truthTable}}}
\examples{
\dontrun{
tt <- fs_tt(Lipset_fs,"Survived.FZ",c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ"),cases="Country",consistency_cutoff=0.7)
eqmcc(tt,"OUT",expl.1 =T,incl.rem=F,show.cases=TRUE,use.letters=FALSE) ## Formula 1 in page 112
eqmcc(tt,"OUT",expl.1 =T,incl.rem=T,show.cases=TRUE,use.letters=FALSE) ## Formula 2 in page 114
eqmcc(tt,"OUT",expl.0 =T,incl.rem=F,show.cases=TRUE,use.letters=FALSE) ## Formula 5 in page 115
eqmcc(tt,"OUT",expl.0 =T,incl.rem=T,show.cases=TRUE,use.letters=FALSE) ## Formula 6 in page 117
}
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
