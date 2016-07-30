\name{suffnec}
\alias{suffnec}
\alias{print.suffnec}
\title{Produces a sufficiency and necessity matrix}
\description{
creates a matrix of the sufficiency and necessity scores for a crisp set
or fuzzy set data frame.
}
\usage{
suffnec(x, use = c("complete", "pairwise"))
%# print method of class suffnec
\method{print}{suffnec}(x, digits = 3, ...)
}
\arguments{
  \item{x}{For coincide, x is a dataframe of crip set or fuzzy set data, which ranges from
  0 to 1.}
  \item{use}{The method of handling missing data. "complete" means
    listwise deletion and "pairwise" means pairwise deletion.}
  \item{digits}{a non-null value for 'digits' specifies the minimum
    number of significant digits to be printed in values.}
  \item{\ldots}{arguments passed to default method of print.}
}
\details{
In the terminology of set theory, if X is the sufficient condition of Y,
then the score of X should be consistently less or equal than that of
Y. Similarly, if X is the necessary condition of Y, then
the sore of X should be consistently greater or equal than that of
Y. The necessary score measures such consitency. The formulas can be found in Ragin(2006:297).

For crip set, sufficient score measures proportation of Y=1 given that X=1. 
The neccessary score measures the proportation of X=1 given that Y=1.
 }
\value{
  A list of two matrixs of consistency score.
  \item{suff}{Sufficiency Scores Matrix,measuring the consistency score
    of 'X is sufficient condition of Y'.}
  \item{nec}{Necessity Scores Matrix, measuring the consistency score of
  'X is necessary condition of Y'.}
}
\references{
Ragin, Charles C. 2006. "Set Relations in Social Research:Evaluating
Their Consistency and Coverage." Political Analysis 14 (3) : 291-310
}

 \author{Ronggui HUANG}
 \seealso{ \code{\link[QCA3]{coincid}}}
 \examples{
suffnec(CoD[,c('cod',"econdev","ethlihom","nocom")]) 
}
