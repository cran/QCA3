\name{fsplot}
\alias{fsplot}
\title{
  Fuzzy set plot.
}
\description{
  Fuzzy set plot facilitating assessment of a causal recipe.
}
\usage{
fsplot(formula, data, main = "fuzzy set plot", xlab = NULL, ylab = NULL,...)
}
\arguments{
  \item{formula}{
    a formula like out ~ condition1 + condition2 + conditions3, with out
    being the outcome and condition1, condition2 and conditions3 are
    conditions in a causal recipe.
  }
  \item{data}{
    a data frame.
  }
  \item{main}{
    a title of the plot.
  }
  \item{xlab}{
    a title for the x axis
  }
  \item{ylab}{
    a title for the y axis
  }
  \item{\dots}{ other arguments passed to \code{plot}.}
}
\details{
  This function is useful for exploratory analysis. If a recipe is a
  causal condition of the outcome, then the fuzzy set score should be
  consistently less or equal than that of outcome. That is, if the
  recipe does not occur, the outcome may occur; however, when the recipe
  occurs, the outcome MUST occurs. It suggests that if a recipe can be
  interpreted as causal condition, the points scatter above the diagonal
  by and large. The set-theoretic consistency measures this
  consitency. High consistency score endors the causal interpretation.

  In addition, high consistency permits interpretation of coverage
  score. This score indicates the proportion of membership in outcome
  accounted for by the recipe. It should noted that it is harzardous to
  interpret the coverage measure when consistency score is low.
}
\references{
  Ragin. C. 2008. "Configurational Thinking." Chapter 6 of "Redesigning
  Social Inquiry." The University of Chicago Press.  
  }
\author{
Ronggui HUANG
}
\seealso{
See Also as \code{\link{coverage}}, \code{\link{consistency}}
}
\examples{
fsplot(cod~econdev,data=CoD)
fsplot(cod~eduhi,data=CoD) ## low consistency, cautious when interpret coverage
fsplot(cod~econdev+eduhi+ethlihom+close,data=CoD)
}