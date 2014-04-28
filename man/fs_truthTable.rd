\name{fs_truthTable}
\alias{fs_truthTable}
\alias{sort.fs_truthTable}
\title{Construction of truthTable from fuzzy set score}
\description{
  Constructing a truthTable from fuzzy set score.
}
\usage{
fs_truthTable(mydata, outcome, conditions, ncases_cutoff = 1,
              consistency_cutoff = 0.8, prop_cutoff=1, 
              show.cases =TRUE,quiet = FALSE,
              cases = NULL, complete=FALSE,...)

\method{sort}{fs_truthTable}(x, decreasing = TRUE, criterion = "Consistency", ...)
}
\arguments{
  \item{mydata}{ A fuzzy set score dataset. All the scores must range
    from 0 to 1.}
  \item{outcome}{ character, the name of the outcome variable in the dataset.}
  \item{conditions}{ character vetor, the name of the conditions from the dataset.}
  \item{ncases_cutoff}{When number of case is less than cutoff, it will
    be regarded as dontcare configuration.}
  \item{consistency_cutoff}{Cutoff point of consistenty score, cases with
    consistency score greater than cutoff point are regarded as OUT=1. }
  \item{prop_cutoff}{proportion of consistent/inconsistent cases.}
  \item{show.cases}{show the rownames from the original dataset for each
    combination of conditions.}
  \item{quiet}{Not used currently.}
  \item{cases}{ character, variable of case names. If it is NULL and
    show.cases is TRUE, name of cases are derived from row names of
    dataset.}
  \item{complete}{show logical remainders when TRUE.}
  \item{\dots}{Not used currently.}

  \item{x}{an fs_truthTable object.}
  \item{decreasing}{same as that of sort.}
  \item{criterion}{a name from the truthTable, sort the fs_truthTable
  according to this variable.}
}
\details{
  There are serveral pillars which make it possible to construct a crip
  truthTable summarizing the raw data. There is a correspondance between
  vector space corners and truthTable rows. Thus, it is possible to get
  the number of cases with strong membership in each corner (usually
  greater then 0.5), and the consistency of the empirical evidence for
  each corner. By specifying the frequency thresholds for fuzzy-set
  assessments (the \code{ncases_cutoff} argument), and assessing the
  consistency of fuzzy-set subset relations (the \code{consistency_cutoff}
  argument), we can finally construct a truthTable.

  This function can also be used for crip set QCA as fsQCA software does,
  though cs_truthTable is specific to csQCA.

  There is a sort method and a consistGaps method for the fs_truthTable object.
}
\value{
  An object of class "fs_truthTable" and "truthTable".
}
\references{
  Ragin. Charles. 2009. Qualitative Comparative Analyais Using Fuzzy
  Sets (fsQCA). In Configuraional comparative Methods: qualitative
  comparative analysis (QCA) and related techniques. ed by Benoit
  RiHoux and Charles Ragin. Sage. This chapter can be downloaded from
  \url{http://www.u.arizona.edu/~cragin/fsQCA/software.shtml}.

  Rubinson, C. (2013), 'Contradictions in fsQCA', Qual Quant 47(5), 2847-2867.
}
\author{Ronggui HUANG}
\seealso{\code{\link{reduce}}, \code{\link{cs_truthTable}} and  \code{\link{fs_truthTable}}}
\examples{
fs_truthTable(Lipset_fs,"Survived.FZ",
             c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ"),
             cases="Country",consistency_cutoff=0.7)

fst <- fs_truthTable(Lipset_fs,"Survived.FZ",
             c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ"),
             cases="Country",consistency_cutoff=0.7,complete=TRUE)

sort(sort(fst),criterion="OUT")

consistGap(fst)
}