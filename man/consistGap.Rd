\name{consistGap}
\alias{consistGap}
\title{
Consistency gaps between configurations
}
\description{
Show the consistency gaps between configurations
}
\usage{
consistGap(x)
}
\arguments{
  \item{x}{
 A fuzzy set truthTable.}
}
\details{
Sometimes, QCA users should not mechanically use a pre-determined consistency cut-off. 
One possible way to determine such a cut-off is to examine the biggest consistency gap 
between configurations (Schneider & Wagemann,2010:407) and place a cutting point between 
the two configuations.

The truthTable is sorted according to consistency and outcome first before calculating 
the consistency gap.
}
\value{
A data frame with 6 columns, with ConsistGap the most important variable.}
\references{
Schneider, C. & Wagemann, C. 2010. Standards of Good Practice in Qualitative Comparative 
Analysis (QCA) and Fuzzy-Sets. Comparative Sociology, 9 (3): 397-418.
}
\author{
Ronggui HUANG}
\seealso{
\link{setOUT}, \link{fs_truthTable}
}
\examples{
xs<-c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ","Stable.FZ")

fst <- fs_truthTable(Lipset_fs,"Survived.FZ", xs,
        cases="Country",consistency_cutoff=0.7)

consistGap(fst)
}