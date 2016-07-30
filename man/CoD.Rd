\name{CoD}
\alias{CoD}
\docType{data}
\title{
The consolidation of Democracy.
}
\description{
Data used in "The consolidation of democracy: comparing Europe and Latin
America", Schneider (2009).
}
\usage{data(CoD)}
\format{
  A data frame with 32 observations on the following 20 variables.
  \describe{
    \item{\code{country}}{a factor of countray names}
    \item{\code{label}}{a factor of country labels}
    \item{\code{CODPERCENTAGE}}{a numeric vector}
    \item{\code{cod}}{a numeric vector}
    \item{\code{GDP}}{a numeric vector}
    \item{\code{econdev}}{a numeric vector}
    \item{\code{EDUCATION}}{a numeric vector}
    \item{\code{eduhi}}{a numeric vector}
    \item{\code{ETHNOLING}}{a numeric vector}
    \item{\code{ethlihom}}{a numeric vector}
    \item{\code{DIST}}{a numeric vector}
    \item{\code{close}}{a numeric vector}
    \item{\code{DEMYEARS}}{a numeric vector}
    \item{\code{demex}}{a numeric vector}
    \item{\code{COMYEARS}}{a numeric vector}
    \item{\code{nocom}}{a numeric vector}
    \item{\code{parlia}}{a numeric vector}
    \item{\code{PARTYNUMBER}}{a numeric vector}
    \item{\code{efpahi}}{a numeric vector}
    \item{\code{decent}}{a numeric vector}
  }
}
%\details{}
\source{
This data set is kindly provided by Prof Schneider.
}
\references{
Schneider, 2009, "The consolidation of democracy: comparing Europe and Latin America". Routledge.
}
\examples{
data(CoD)
# Examples from Chapter 7
### table 6.2 in p77
cond <- c("econdev","eduhi","ethlihom","close","demex","nocom","parlia","efpahi","decent")
round(apply(CoD[, cond], 2, FUN=function(x) consistency(x,CoD$cod,'great')),3)
round(apply(CoD[, cond], 2, FUN=function(x) consistency(1-x,CoD$cod,'great')),3)
round(apply(CoD[, cond], 2, FUN=function(x) consistency(x,1-CoD$cod,'great')),3)
round(apply(CoD[, cond], 2, FUN=function(x)
    
consistency(1-x,1-CoD$cod,'great')),3)

round(apply(CoD[,c("econdev","demex")], 2, FUN=function(x) coverage(1-x,1-CoD$cod,'great')),3)
round(consistency(apply(1-CoD[,c("econdev","demex")],1,min),1-CoD$cod,'great'),3)
round(coverage(apply(1-CoD[,c("econdev","demex")],1,min),1-CoD$cod,'great'),3)

## QCA of remote conditions
reduce(CoD,'cod',c('econdev','ethlihom','close','demex','nocom'),'positive','include',
       case='label',pre='fs_truthTable',keep=FALSE)

## consitency of remote conditions as sufficent condition (table 6.4 in p82)
round(apply(CoD[,c("econdev","ethlihom","nocom")], 2, 
           FUN=function(x) consistency(x,CoD$cod,'less')),3)
           
round(apply(CoD[,c("econdev","ethlihom","nocom")], 2, 
          FUN=function(x) coverage(x,CoD$cod,'less')),3)
          
round(consistency(apply(CoD[,c("econdev","ethlihom","nocom")],1,max),CoD$cod),3)
round(coverage(apply(CoD[,c("econdev","ethlihom","nocom")],1,max),CoD$cod),3)

## table 6.4 in page 82
fs=fs_truthTable(CoD,'cod',c('econdev','eduhi','ethlihom','close','demex','nocom'),case='label')
fs <- sort(fs,,"Consistency")
tb6.4 <- reduce(fs,remainder='include')

## table 6.6 in page 86
reduce(CoD,'cod', cond=c('parlia','efpahi','decent'),'positive','include',
       case='label',pre='fs_truthTable')
       
reduce(CoD,'cod',c('econdev','ethlihom','nocom','parlia','efpahi','decent'),'positive','include',
       case='label',pre='fs_truthTable')

## table 6.7 in pp87-88
fs_truthTable(CoD,'cod',c('econdev','ethlihom','nocom','parlia','efpahi','decent'),case='label')
}
\keyword{datasets}
