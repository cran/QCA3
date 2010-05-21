\name{cs_truthTable}
\alias{cs_truthTable}
\alias{csTQCA_truthTable}
\title{Construct a truthTable for csQCA}
\description{
  Construct a truthTable for csQCA. Both deterministic and
  probabilistic methods of determining configurations of positive,
  negative and contraditory outcome are implemented.
  
  This function can be used for crip set TQCA as well. See
  \code{\link{CarenPanofsky}} for example. It needs manual construction
  of indicator conditions of temporal order.
}
\usage{
cs_truthTable(mydata, outcome, conditions,method = c("deterministic", 
    "probabilistic"), weight = NULL, show.cases = TRUE, cases = NULL, 
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
  The value of all the conditions should start from 0. In
  cript set QCA, it is always be 0 or 1. Value -9 in conditions means "don't care"
  (though "don't care" in outcome is denoted by "-9").
  
  Symbols of the outcome. '1' is postive configuration, '0' is negative
  configuration, 'C' is contraditory configuration, "?" is unobserved
  configuration and '-9' is don't care configuration. When show.case is
  TRUE and a configuration is 'C', then the name of case with negative
  outcome is highlighted by [].
  
  'cutoff1' and 'cutoff0' are only meaningful for'deterministic'
  method. They represent cutting point of positive case and negative
  case. When a configuration has positive case only and the number of
  cases is equal or greater than the cutting point, then it is regared as
  positive outcome, otherwise as don't care outcome. Similarly, When a
  configuration has negative case only and the number of cases is equal or
  greater than the cutting point, then it is regared as negative
  outcome, otherwise as don't care outcome. If a configuration has both
  positive case and nagetive case, the number of each type of cases will
  be compared with the corresponding cutting point. If only one type of 
  case have enough case (number of cases greater than cutting point),
  that configuration is regarded as that type. If both types have enough
  case, it is contraditory configuration. If neither type has enough
  case, it is don't care configuration.

  The caculation of cutting point: if it is equal or greater than 1, the
  cutting point is the value of cutoff1 and cutoff0. If it is between 0
  and 1, then the cutting point is the cutoff1/cutoff0 multiplied by the
  total number of case.

  'benchmark' and 'conf.level' are only meaningful for 'probabilistic'
  method. When the method is 'probabilistic', a statistical test will
  conducted to test if the proportion of case for a configuration is
  greater then a benchmark. If the proportion of cases with outcome '1' is
  greater than benchmark, then the it is a configuratin with outcome
  '1'. If the proportion of case with outcome '0' is greater than
  benchmark, then the configuration with outcome of '0'. If neither
  proportion fits the criterion, then it is don't care
  configuration. There is no contraditory congfiguration in this method,
  as it is designed to handle with contraditory configurations.

  rownames of a truthTable is grouping index (not important for end-users).
}
\value{
  An object of class "truthTable" and "cs_trutbTable". A list with 5
  components:
  \item{truthTable}{a data frame presenting a truthTable.}
  \item{outcome}{The name of outcome variable. length-1 character.}
  \item{conditions}{conditions. A character vector.}
  \item{nlevels}{integer vector specifying number of levels of each condition.}
  \item{call}{the matched call.}
}
\references{
  Ragin, Charles. 2000. Fuzzy-Set Social Science. Pp109-116. University Of Chicago Press.
}
\author{Ronggui HUANG}
\seealso{\code{\link{fs_truthTable}},  \code{\link{mv_truthTable}} and \code{\link{reduce}}}
\examples{
## truthTable for csQCA
cs_truthTable(Lipset_cs,"SURVIVAL", c("GNPCAP", "URBANIZA", "LITERACY",
  "INDLAB", "GOVSTAB"),case="CASEID")
}

