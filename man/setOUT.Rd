\name{setOUT}
\alias{setOUT}
\title{
Reset the outcome of configurations in a truthTable
}
\description{
To manually change the outcome of configurations in a truthTable.
}
\usage{
setOUT(x, rownames, value)
}
\arguments{
  \item{x}{a truthTable.}
  \item{rownames}{character vector, rownames of the configurations to be changed.}
  \item{value}{the new vaule,one of '1', '0', '-9' or '?'.}
}
\details{
Sometimes, researcher want to override the outcome for theoretical reasons. 
This function provides a shortcut so the researchers do not have to mess with the internal structure of a truthTable.
This is inspired by fsQCA (2.5) which require the researcher to decide the outcome manually or semi-automatically.

Should note that "freq1" and "freq0" are not affected, so QCA methods reliant on these two pieces of information might not suitable any more.
}
\value{
A new truthTable with a new outcome.
}
\seealso{
\link{consistGap}
}
\author{Ronggui HUANG}
\examples{
fst <- fs_truthTable(Lipset_fs, "Survived.FZ",
                c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ","Stable.FZ"),
                cases="Country",consistency_cutoff=0.7)
setOUT(fst,"240","0")
}