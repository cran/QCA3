\name{mv_truthTable}
\alias{mv_truthTable}
\alias{sort.mv_truthTable}
\title{Construct a truthTable for mvQCA}
\description{
  Construct a truthTable for mvQCA. Both deterministic and
  probabilistic methods of determining configurations of positive,
  negative and contraditory outcome are implemented.
}
\usage{
mv_truthTable(mydata, outcome, conditions, method = c("deterministic", 
    "probabilistic"), weight = NULL, complete=FALSE, show.cases = TRUE, cases = NULL, 
    cutoff1 = 1, cutoff0 = 1, benchmark = 0.65, conf.level = 0.95, 
    missing = c("missing", "dontcare", "positive", "negative"))
}
\arguments{
  \item{mydata}{data frame of the raw data.}
  \item{outcome}{character, the name of the outcome variable in mydata.}
  \item{conditions}{character vector, the name of the conditions from mydata.}
  \item{method}{character, specifying the method of determining the
    outcome of a configuration.}
  \item{weight}{character, name of a variable specifying the weights.}
  \item{complete}{logical, when it is TRUE the result includes
    configurations without empirical cases.}
  \item{show.cases}{logical, when TRUE the result shows case names.}
  \item{cases}{character, variable specifying the case names. When it is
    NUll, then use row names of mydata as case names.}
  \item{cutoff1}{length one numeric vector.}
  \item{cutoff0}{length one numeric vector. }
  \item{benchmark}{Benchmark for statistical test. Must equal or greater
    than 0.5.}
  \item{conf.level}{confident level of statistical test.}
  \item{missing}{method to handle missing data.}
}
\details{
  The value of all the conditions should start from 0. For example, for
  mvQCA with 3 levels,it should be 0,1 or 2. Value -9 in conditions
  means "don't care" (though "don't care" in outcome is denoted by "-9").

  There is a sort method method for the truthTable object.
  
  Refer to \code{\link{cs_truthTable}} for for more details.
}
\value{
  An object of class "mv_truthTable" and "truthTable".
}
\references{
  Ragin, Charles. 2000. Fuzzy-Set Social Science. Pp109-116. University Of Chicago Press.
}
\author{Ronggui HUANG}
\seealso{\code{\link{fs_truthTable}} and \code{\link{reduce}}}
\examples{
## truthTable for mvQCA.
mv_truthTable(Lipset_mv,"SURVIVAL", c("GNPCAP", "URBANIZA", "LITERACY",
              "INDLAB"),case="CASEID")

mvt <- mv_truthTable(Lipset_mv,"SURVIVAL", c("GNPCAP", "URBANIZA", "LITERACY",
              "INDLAB"),case="CASEID")

sort(mvt)
}
