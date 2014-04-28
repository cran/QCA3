\name{QCA.methods}
\alias{print.QCA}
\alias{summary.QCA}
%\alias{plot.QCA}
\alias{[.QCA}
\alias{update.QCA}
\title{Methods for "QCA" an object}
\description{
   Various methods for object from \code{\link{reduce}}
}
\usage{
%\method{plot}{QCA}(x,...)
%
\method{print}{QCA}(x, traditional = TRUE, show.truthTable = FALSE, ...)

\method{summary}{QCA}(object, traditional = TRUE, show.case = TRUE, ...)

\method{[}{QCA}(object, which)

\method{update}{QCA}(object,...,evaluate = TRUE)
}
\arguments{
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
  \item{\dots}{ For \code{print.QCA} and \code{summary.QCA}, currently not
    used. For \code{update}, additional arguments to the call, or
arguments with changed values. Use 'name=NULL' to remove the argument
'name'.}
\item{evaluate}{when TRUE, return the evaluated result which is an
object of QCA class. Otherwise, it returns the call.}
}
\details{
  The traditional way uses upper-case letters representing 1 and and
  lower-case letters reprensenting 0. The Tosmana-style uses
  \code{condition{value}} to represent the prime implicants.

  The summary method calculates the number of cases covered by each implicant,
  the percentage of explained cases, the number of cases (percentage of explained cases)
  covered by multiple
  implications.It also shows the case names covered by each implicant, with
  (n) indicating the number of implicants coving the case(s).
}
\value{
 print method does not return any value.

 summary method returns an object of class "summary.QCA".

 The index method returns an object of class "QCA".

 update method returns a new "QCA" object if evaluate is TRUE, the call
if FALSE.
}

\author{ Ronggui HUANG}

\seealso{
\code{\link{reduce}} and \code{\link{constrReduce}}
}

