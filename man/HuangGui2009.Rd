\name{HuangGui2009}
\Rdversion{1.1}
\alias{HuangGui2009}
\docType{data}
\title{
data about on-line forum and homeowners' collective resistance
}
\description{Dataset used in Huang and Gui's (2009) study of Internet and homeowners' collective resistance.
}
\usage{data(HuangGui2009)}
\format{
  A data frame with 15 observations on the following 6 variables.
  \describe{
    \item{\code{network}}{stock of social ties.}
    \item{\code{big}}{scale of neighborhood.}
    \item{\code{grievance}}{Existence of grevance.}
    \item{\code{ha}}{Existence of Homeowners' Association.}
    \item{\code{BBS}}{Existence of homeowners' onl-line forum.}
    \item{\code{outcome}}{outcome, existence of collective resistance.}
    \item{\code{code}}{id of neighborhood.}
  }
}
%\details{}
\source{
HUANG Ronggui, GUI Yong, 2009, Internet and homeowners' collective
resistance: a comparative qualitative analysis in
Shanghai". Sociological Research.(5):29-56.(In Chinese)
}
\examples{
data(HuangGui2009)
reduce(HuangGui2009,"outcome",c("network","big","grievance","ha","BBS"),"pos",case="code")
ans1 <- reduce(HuangGui2009,"outcome",c("network","big","grievance","ha","BBS"),"pos",case="code",rem="include")
SA(ans1)
reduce(HuangGui2009,"outcome",c("network","big","grievance","ha","BBS"),"neg",case="code")
ans2 <- reduce(HuangGui2009,"outcome",c("network","big","grievance","ha","BBS"),"neg",case="code",rem="incl")
CSA(ans1,ans2[3])
SA(ans2[3])
}
\keyword{datasets}
