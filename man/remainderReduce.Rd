\name{remainderReduce}
\alias{remainderReduce}
\alias{print.remainders}
\title{
Remainder configurations.
}
\description{
Minimizing the remainder configurations.
}
\usage{
remainderReduce(x, keepTruthTable = TRUE)
}
\arguments{
  \item{x}{
A truthTable object, from \link{cs_truthTable}, \link{fs_truthTable} or \link{mv_truthTable}
}
  \item{keepTruthTable}{ logical, when TRUE the returned object keeps the truthTable
}
}
%\details{}
\value{
An object of class c("remainders","QCA"), the components are the same as that from \link{reduce}.
}
%%\references{}
\author{Ronggui HUANG}
\examples{
fst <- fs_truthTable(Lipset_fs,"Survived.FZ",
       c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ"),
       cases="Country",consistency_cutoff=0.7,complete=TRUE)
## complete must be TRUE
remainderReduce(fst)
}