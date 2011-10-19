\name{CoD}
\alias{CoD}
\docType{data}
\title{
The consolidation of Democracy.
}
\description{
Data used in "The consolidation of democracy: comparing Europe and Latin America", Schneider (2009).
}
\usage{data(CoD)}
\format{
  A data frame with 32 observations on the following 20 variables.
  \describe{
    \item{\code{country}}{a factor with levels \code{Albania} \code{Argentina} \code{Belarus} \code{Bolivia} \code{Brazil} \code{Bulgaria} \code{Chile} \code{Czech_R} \code{Ecuador} \code{Estonia} \code{Georgia} \code{Greece} \code{Guatemala} \code{Honduras} \code{Hungary} \code{Latvia} \code{Lithuania} \code{Mexico} \code{Mongolia} \code{Nicaragua} \code{Paraguay} \code{Peru} \code{Poland} \code{Portugal} \code{Romania} \code{Russia} \code{Slovakia} \code{Slovenia} \code{Spain} \code{Turkey} \code{Ukraine} \code{Uruguay}}
    \item{\code{label}}{a factor with levels \code{AL} \code{AR} \code{BE} \code{BO} \code{BR} \code{BU} \code{CH} \code{CR} \code{EC} \code{EST} \code{GE} \code{GR} \code{GUA} \code{HO} \code{HU} \code{LAT} \code{LIT} \code{MO} \code{MX} \code{NI} \code{PA} \code{PE} \code{PL} \code{PO} \code{RO} \code{RU} \code{SK} \code{SL} \code{SP} \code{TU} \code{UK} \code{UR}}
    \item{\code{CODPERCENTAGE}}{a factor with levels \code{14.58} \code{26.14} \code{26.39} \code{31.48} \code{33.58} \code{39.46} \code{43.45} \code{48.96} \code{49.34} \code{54.17} \code{54.66} \code{62.5} \code{65.42} \code{65.91} \code{66.67} \code{68.33} \code{69.17} \code{70.14} \code{70.45} \code{71.25} \code{78.55} \code{78.85} \code{79.92} \code{82.29} \code{83.58} \code{missing}}
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
round(apply(CoD[,c("econdev","eduhi","ethlihom","close","demex","nocom","parlia","efpahi","decent")],2,FUN=function(x) consistency(x,CoD$cod,'great')),3)
round(apply(CoD[,c("econdev","eduhi","ethlihom","close","demex","nocom","parlia","efpahi","decent")],2,FUN=function(x) consistency(1-x,CoD$cod,'great')),3)
round(apply(CoD[,c("econdev","eduhi","ethlihom","close","demex","nocom","parlia","efpahi","decent")],2,FUN=function(x) consistency(x,1-CoD$cod,'great')),3)
round(apply(CoD[,c("econdev","eduhi","ethlihom","close","demex","nocom","parlia","efpahi","decent")],2,FUN=function(x) consistency(1-x,1-CoD$cod,'great')),3)
round(apply(CoD[,c("econdev","demex")],2,FUN=function(x) coverage(1-x,1-CoD$cod,'great')),3)
round(consistency(apply(1-CoD[,c("econdev","demex")],1,min),1-CoD$cod,'great'),3)
round(coverage(apply(1-CoD[,c("econdev","demex")],1,min),1-CoD$cod,'great'),3)

## QCA of remote conditions
reduce(CoD,'cod',c('econdev','ethlihom','close','demex','nocom'),'positive','include',case='label',pre='fs_truthTable',keep=FALSE)

## consitency of remote conditions as sufficent condition (table 6.4 in p82)
round(apply(CoD[,c("econdev","ethlihom","nocom")],2,FUN=function(x) consistency(x,CoD$cod,'less')),3)
round(apply(CoD[,c("econdev","ethlihom","nocom")],2,FUN=function(x) coverage(x,CoD$cod,'less')),3)
round(consistency(apply(CoD[,c("econdev","ethlihom","nocom")],1,max),CoD$cod),3)
round(coverage(apply(CoD[,c("econdev","ethlihom","nocom")],1,max),CoD$cod),3)

## table 6.4 in page 82
fs=fs_truthTable(CoD,'cod',c('econdev','eduhi','ethlihom','close','demex','nocom'),case='label')
fs <- sort(fs,,"Consistency")
tb6.4 <- reduce(fs,remainder='include')

## table 6.6 in page 86
reduce(CoD,'cod',c('parlia','efpahi','decent'),'positive','include',case='label',pre='fs_truthTable')
reduce(CoD,'cod',c('econdev','ethlihom','nocom','parlia','efpahi','decent'),'positive','include',case='label',pre='fs_truthTable')
## table 6.7 in pp87-88
fs_truthTable(CoD,'cod',c('econdev','ethlihom','nocom','parlia','efpahi','decent'),case='label')

}
\keyword{datasets}
