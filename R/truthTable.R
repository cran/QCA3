## This file (R/truthTable.R) is part of QCA3 package
## copyright: HUANG Ronggui 2008-2010

lowerLimite <- function(x, n, conf.level=0.95) {
  ## If lowerLimite() > benchmark, then it the result supportes H1.
  ## see Ragin (2001:109-115).
  ans <- numeric(length(x))
  idx <- which(x!=0)
  ans[idx] <- qbeta(1 - conf.level, x[idx], n[idx] - x[idx] + 1)
  ans
}

cs_truthTable <- function(mydata, outcome, conditions,
                          method = c("deterministic","probabilistic"),
                          weight=NULL,
                          show.cases = TRUE, cases=NULL,
                          cutoff1 = 1, cutoff0 = 1, benchmark=0.65, conf.level = 0.95,
                          missing=c('missing','dontcare','positive','negative')
                          )
{
    ## mydata is a data frame.
    ## outcome and conditions are character vectors.
    nlevels = rep(2,length(conditions))
    if (outcome==""||conditions =="") stop("You must specific outcome and conditions first.")
    if (length(conditions)<2) stop("The number of conditions must greater than 1.")
    reserved <- c("NCase","freq1","freq0","OUT","Cases")
    if (any(outcome %in% reserved)) stop("Some names of condition are reserved for truthTable.")
    mydata <- mydata[,c(outcome,conditions,weight,cases)]
    missing <- match.arg(missing)
    if (missing=="missing")  mydata <- na.exclude(mydata) # eliminate missing data
    if (missing=='dontcare') mydata[is.na(mydata)] <- -9
    if (missing=='positive') mydata[is.na(mydata)] <- 1
    if (missing=='negative') mydata[is.na(mydata)] <- 0
    ## take care of missing data
    fulldata <- mydata[,c(outcome,conditions)]
    outcomeData <- fulldata[,outcome]
    if (any(!outcomeData %in% c(0,1))) stop("outcome value must in [0,1].")
    conditionsData <- fulldata[,conditions]
    colmax <- sapply(conditionsData,max,na.rm=T)
    if (any(colmax+1 > nlevels)) {
        warning("It seems multi-value QCA, use 'mv_truthTable' instead of 'cs_truthTable'.")
        nlevels <- colmax + 1
    }
    if (!is.null(weight)) weight <- mydata[[weight]] else weight <- rep(1, nrow(mydata))
    method <- match.arg(method)
    ## getId <- function(implicant,nlevels){
    ##     ## id of combinations
    ##     IDX <- cumprod(nlevels)/nlevels
    ##     ans <- sum(implicant*IDX)+1
    ##     ans
    ## }
    ## rowid <- apply(conditionsData, 1, getId, nlevels=nlevels)
    rowid <- apply(conditionsData, 1, implicant2Id, nlevels=nlevels)
    ## use id of grouping rather than combination to handle dontcare case
    N_total <- sum(weight,na.rm=TRUE) ## total number of case taking freq weight into consideration
    Positive <- tapply(outcomeData,rowid,FUN=function(each) all(each==1)) ## index of configuration with positive outcome
    ## with aid of rowid, we aggregate those with common rowid into one group.
    Pid <- names(Positive)[Positive] ## rownames of configuration with positive outcome
    Negative <- tapply(outcomeData,rowid,FUN=function(each) all(each==0))
    Nid <- names(Negative)[Negative]
    Contradictory <- tapply(outcomeData,rowid,FUN=function(each) {
        c1 <- (!all(each==0)) && (!all(each==1))
        c1})
    Cid <- names(Negative)[Contradictory] ## all.equal(names(Positive),names(Negative))
    WhichUnique <- match(sort(unique(rowid)),rowid) ## pay attention to the use of match()
    allExpress <- conditionsData[WhichUnique,]
    rownames(allExpress) <- as.character(sort(unique(rowid)))
    ## NCase
    allExpress$NCase <- 0
    Ncase <- tapply(weight,rowid,sum)
    allExpress$NCase[match(names(Ncase),rownames(allExpress))] <- Ncase
    ## freq0 and freq1
    ## allExpress$freq0 <- allExpress$freq1 <- "-9"
    allExpress$freq0 <- allExpress$freq1 <- 0
    Ncase1 <- by(cbind(weight,outcomeData),rowid,FUN=function(idx) sum(idx[,1][idx[,2]==1]))
    allExpress$freq1[match(names(Ncase1),rownames(allExpress))] <- Ncase1
    Ncase0 <- by(cbind(weight,outcomeData),rowid,FUN=function(idx) sum(idx[,1][idx[,2]==0]))
    allExpress$freq0[match(names(Ncase0),rownames(allExpress))] <- Ncase0
    ## out status
    allExpress$OUT <- "?"
    if (method=="deterministic"){
        cutoff1 <- ifelse(cutoff1<1,cutoff1*N_total,cutoff1)
        cutoff0 <- ifelse(cutoff0<1,cutoff0*N_total,cutoff0)
        pidx <- intersect(match(Pid,rownames(allExpress)), which(allExpress$freq1 >= cutoff1))
        allExpress$OUT[pidx] <- "1"
        nidx <- intersect(match(Nid,rownames(allExpress)),which(allExpress$freq0 >= cutoff0))
        allExpress$OUT[nidx] <- "0"
        cidx1 <- intersect(match(Cid,rownames(allExpress)),which(allExpress$freq1 >= cutoff1))
        cidx0 <- intersect(match(Cid,rownames(allExpress)),which(allExpress$freq0 >= cutoff0))
        cidx <- intersect(cidx1, cidx0)
        allExpress$OUT[cidx]<-"C"
        Dontcare1 <- intersect(match(Cid,rownames(allExpress)),which(allExpress$freq1 < cutoff1))
        Dontcare0 <- intersect(match(Cid,rownames(allExpress)),which(allExpress$freq0 < cutoff0))
        Dontcareid <- intersect(Dontcare1, Dontcare0)
        allExpress$OUT[Dontcareid]<-"-9"
        allExpress$OUT[intersect(cidx1,Dontcare0)]<-"1"
        allExpress$OUT[intersect(cidx0,Dontcare1)]<-"0"
        ## Dontcareid <- as.character(setdiff(rowid,rownames(allExpress)[c(pidx,nidx,cidx)]))
        ## with Ncases less then cutoff point.
        ## allExpress$OUT[match(Dontcareid,rownames(allExpress))] <- "-9"
    }
    if (method=="probabilistic"){
        limit1 <- lowerLimite(allExpress$freq1,allExpress$NCase,conf.level)
        limit0 <- lowerLimite(allExpress$freq0,allExpress$NCase,conf.level)
        pidx <- intersect(which(limit1 >=benchmark),match(c(Pid,Cid),rownames(allExpress)))
        nidx <- intersect(which(limit0 >=benchmark),match(c(Nid,Cid),rownames(allExpress)))
        Dontcareid <- setdiff(match(c(Nid,Cid,Pid),rownames(allExpress)),c(pidx,nidx))
        allExpress$OUT[pidx] <- "1"
        allExpress$OUT[nidx] <- "0"
        allExpress$OUT[Dontcareid] <- "-9"
        ## no contradictory cases when using probabilistic method???
    }
    ## show.cases
    if (show.cases){
        if (is.null(cases)) casesNames <- rownames(mydata) else casesNames <- mydata[,cases]
        casesNames <- gsub(",","_",casesNames)
        casesNames[outcomeData==0] <- paste("[",casesNames[outcomeData==0],"]",sep="") ## mark the negative cases
        casesNames <- tapply(casesNames,rowid,FUN=function(each) paste(each,sep="",collapse=", "))
        allExpress$Cases <- ""
        allExpress$Cases[match(names(casesNames),rownames(allExpress))] <- casesNames
        allExpress$Cases[allExpress$OUT!="C"] <- gsub("\\[|\\]","",allExpress$Cases[allExpress$OUT!="C"]) ## mark contr case
    }
    allExpress
    ans <- list(truthTable=allExpress,outcome=outcome,conditions=conditions,nlevels=nlevels,call=match.call())
    class(ans) <- c("truthTable","cs_truthTable")
    ans
}



mv_truthTable <- function(mydata, outcome, conditions,
                          method = c("deterministic","probabilistic"),
                          weight=NULL,
                          show.cases = TRUE, cases=NULL,
                          cutoff1 = 1, cutoff0 = 1, benchmark=0.65, conf.level = 0.95,
                          missing=c('missing','dontcare','positive','negative')
                          )
{
    ## mydata is a data frame.
    ## outcome and conditions are character vectors.
    nlevels <- sapply(mydata[,conditions], function(x) max(x,na.rm = T)+1)
    if (outcome==""||conditions =="") stop("You must specific outcome and conditions first.")
    if (length(conditions)<2) stop("The number of conditions must greater than 1.")
    reserved <- c("NCase","freq1","freq0","OUT","Cases")
    if (any(outcome %in% reserved)) stop("Some names of condition are reserved fro truthTable.")
    mydata <- mydata[,c(outcome,conditions,weight,cases)]
    missing <- match.arg(missing)
    if (missing=="missing")  mydata <- na.exclude(mydata) # eliminate missing data
    if (missing=='dontcare') mydata[is.na(mydata)] <- -9
    if (missing=='positive') mydata[is.na(mydata)] <- 1
    if (missing=='negative') mydata[is.na(mydata)] <- 0
    ## take care of missing data
    fulldata <- mydata[,c(outcome,conditions)]
    outcomeData <- fulldata[,outcome]
    if (any(!outcomeData %in% c(0,1))) stop("outcome value must in [0,1].")
    conditionsData <- fulldata[,conditions]
    colmax <- sapply(conditionsData,max,na.rm=T)
    if (any(colmax+1 > nlevels)) {
        warning(sprintf("Mismatch of values of conditions and 'nlevels' argument. \n Replace it with possible value c(%s)",paste(colmax+1,collapse=",")))
        nlevels <- colmax + 1
    }
    if (!is.null(weight)) weight <- mydata[[weight]] else weight <- rep(1, nrow(mydata))
    method <- match.arg(method)
    rowid <- apply(conditionsData, 1, implicant2Id, nlevels=nlevels)
    ## use id of grouping rather than combination to handle dontcare case
    N_total <- sum(weight,na.rm=TRUE) ## total number of case taking freq weight into consideration
    Positive <- tapply(outcomeData,rowid,FUN=function(each) all(each==1)) ## index of configuration with positive outcome
    ## with aid of rowid, we aggregate those with common rowid into one group.
    Pid <- names(Positive)[Positive] ## rownames of configuration with positive outcome
    Negative <- tapply(outcomeData,rowid,FUN=function(each) all(each==0))
    Nid <- names(Negative)[Negative]
    Contradictory <- tapply(outcomeData,rowid,FUN=function(each) {
        c1 <- (!all(each==0)) && (!all(each==1))
        c1})
    Cid <- names(Negative)[Contradictory] ## all.equal(names(Positive),names(Negative))
    WhichUnique <- match(sort(unique(rowid)),rowid) ## pay attention to the use of match()
    allExpress <- conditionsData[WhichUnique,]
    rownames(allExpress) <- as.character(sort(unique(rowid)))
    ## NCase
    allExpress$NCase <- 0
    Ncase <- tapply(weight,rowid,sum)
    allExpress$NCase[match(names(Ncase),rownames(allExpress))] <- Ncase
    ## freq0 and freq1
    allExpress$freq0 <- allExpress$freq1 <- 0
    Ncase1 <- by(cbind(weight,outcomeData),rowid,FUN=function(idx) sum(idx[,1][idx[,2]==1]))
    allExpress$freq1[match(names(Ncase1),rownames(allExpress))] <- Ncase1
    Ncase0 <- by(cbind(weight,outcomeData),rowid,FUN=function(idx) sum(idx[,1][idx[,2]==0]))
    allExpress$freq0[match(names(Ncase0),rownames(allExpress))] <- Ncase0
    ## out status
    allExpress$OUT <- "?"
    if (method=="deterministic"){
        cutoff1 <- ifelse(cutoff1<1,cutoff1*N_total,cutoff1)
        cutoff0 <- ifelse(cutoff0<1,cutoff0*N_total,cutoff0)
        pidx <- intersect(match(Pid,rownames(allExpress)), which(allExpress$freq1 >= cutoff1))
        allExpress$OUT[pidx] <- "1"
        nidx <- intersect(match(Nid,rownames(allExpress)),which(allExpress$freq0 >= cutoff0))
        allExpress$OUT[nidx] <- "0"
        cidx1 <- intersect(match(Cid,rownames(allExpress)),which(allExpress$freq1 >= cutoff1))
        cidx0 <- intersect(match(Cid,rownames(allExpress)),which(allExpress$freq0 >= cutoff0))
        cidx <- intersect(cidx1, cidx0)
        allExpress$OUT[cidx]<-"C"
        Dontcare1 <- intersect(match(Cid,rownames(allExpress)),which(allExpress$freq1 < cutoff1))
        Dontcare0 <- intersect(match(Cid,rownames(allExpress)),which(allExpress$freq0 < cutoff0))
        Dontcareid <- intersect(Dontcare1, Dontcare0)
        allExpress$OUT[Dontcareid]<-"-9"
        allExpress$OUT[intersect(cidx1,Dontcare0)]<-"1"
        allExpress$OUT[intersect(cidx0,Dontcare1)]<-"0"
    }
    if (method=="probabilistic"){
        limit1 <- lowerLimite(allExpress$freq1,allExpress$NCase,conf.level)
        limit0 <- lowerLimite(allExpress$freq0,allExpress$NCase,conf.level)
        pidx <- intersect(which(limit1 >=benchmark),match(c(Pid,Cid),rownames(allExpress)))
        nidx <- intersect(which(limit0 >=benchmark),match(c(Nid,Cid),rownames(allExpress)))
        Dontcareid <- setdiff(match(c(Nid,Cid,Pid),rownames(allExpress)),c(pidx,nidx))
        allExpress$OUT[pidx] <- "1"
        allExpress$OUT[nidx] <- "0"
        allExpress$OUT[Dontcareid] <- "-9"
        ## no contradictory cases when using probabilistic method???
    }
    ## show.cases
    if (show.cases){
        if (is.null(cases)) casesNames <- rownames(mydata) else casesNames <- mydata[,cases]
        casesNames <- gsub(",","_",casesNames)
        casesNames[outcomeData==0] <- paste("[",casesNames[outcomeData==0],"]",sep="") ## mark the negative cases
        casesNames <- tapply(casesNames,rowid,FUN=function(each) paste(each,sep="",collapse=", "))
        allExpress$Cases <- ""
        allExpress$Cases[match(names(casesNames),rownames(allExpress))] <- casesNames
        allExpress$Cases[allExpress$OUT!="C"] <- gsub("\\[|\\]","",allExpress$Cases[allExpress$OUT!="C"]) ## mark contr case
    }
    allExpress
    ans <- list(truthTable=allExpress,outcome=outcome,conditions=conditions,nlevels=nlevels,call=match.call())
    class(ans) <- c("truthTable","mv_truthTable")
    ans
}

fs_truthTable <- function(mydata, outcome, conditions,ncases_cutoff=1,consistency_cutoff=0.8,
                          show.cases = TRUE, quiet = FALSE,cases=NULL,...)
{
    membership_cutoff=0.5
    if (consistency_cutoff>1 || consistency_cutoff<0) stop("consistency_cutoff should be in [0,1].")
    if (consistency_cutoff<0.75) warning("It is suggested that consistency_cutoff be >= 0.75.")
    if (outcome==""||conditions=="") stop("You must specific outcome and conditions first.")
    if (length(conditions)<2) stop("The number of conditions must greater than 1.")
    reserved <- c("NCase","freq1","freq0","OUT","Cases")
    if (any(outcome %in% reserved)) stop("Some names of condition are reserved fro truthTable.")
    mydata <- mydata[,c(outcome,conditions,cases)]
    mydata <- na.exclude(mydata) # eliminate missing data
    fulldata <- mydata[,c(outcome,conditions)]
    if (any(fulldata<0)|| any(fulldata>1)) stop("Fuzzy set score must in [0,1].")
    ncases_cutoff <- ifelse(ncases_cutoff<1,ncases_cutoff*nrow(fulldata),ncases_cutoff)
    allExpress <- eval(parse(text=(sprintf("expand.grid(%s)",paste(conditions,"=1:0",sep="",collapse=",")))))
    conditionsData <- mydata[,conditions]
    ## helper function of getScore
    getScore <- function(index,data){
        Negative <- which(index==0)
        Positive <- which(index==1)
        if (length(Negative)>0 && length(Positive)>0) {
            score <- pmin(apply(1-data[,Negative,drop=FALSE],1,min),apply(data[,Positive,drop=FALSE],1,min))
        } else if (length(Negative)>0 && length(Positive)==0) {
            score <- apply(1-data[,Negative,drop=FALSE],1,min)
        } else if (length(Negative)==0 && length(Positive)>0) {
            score <- apply(data[,Positive,drop=FALSE],1,min)
        }
    }
    ## end of helper function of getScore
    score_mat <- apply(allExpress,1,function(x) getScore(x,data=conditionsData))
    allExpress$NCase<- apply(score_mat,2,function(x) sum(x>membership_cutoff))
    allExpress$Consistency <- apply(score_mat,2,function(x,outcome) {sum(pmin(x,outcome))/sum(x)},outcome=mydata[,outcome])
    allExpress$OUT <- "?"
    allExpress$OUT[allExpress$NCase >= ncases_cutoff & allExpress$Consistency > consistency_cutoff]<-"1"
    allExpress$OUT[allExpress$NCase >= ncases_cutoff & allExpress$Consistency <= consistency_cutoff]<-"0"
    allExpress$OUT[allExpress$NCase < ncases_cutoff & allExpress$NCase >0] <- "-9"
    allExpress$freq0 <- allExpress$freq1 <- 0
    allExpress$freq0[allExpress$OUT=="0"] <- allExpress$NCase[allExpress$OUT=="0"]
    allExpress$freq1[allExpress$OUT=="1"] <- allExpress$NCase[allExpress$OUT=="1"]
    allExpress <- allExpress[,c(seq_len(length(conditions)),(length(conditions)+3):(length(conditions)+5),(length(conditions)+1):(length(conditions)+2))]
    ## reorder alExpress
    if (show.cases){
        if (is.null(cases)) cases <- rownames(mydata) else cases <- mydata[,cases]
        cases <- gsub(",","_",cases)
        allExpress$Cases <- apply(score_mat,2,function(x) paste(cases[which( x > membership_cutoff)],sep="",collapse=","))
    }
    allExpress <- allExpress[allExpress$OUT != "?",,drop=FALSE]
    rownames(allExpress) <- apply(allExpress[,conditions],1, implicant2Id, nlevels=rep(2,length(conditions)))
    ans <- list(truthTable=allExpress,outcome=outcome,conditions=conditions,nlevels=rep(2,length(conditions)),call=match.call())
    class(ans) <- c("truthTable","fs_truthTable")
    ans
}

print.truthTable <- function(x,...){
    x <- unclass(x)
    print(x$truthTable)
}

sort.fs_truthTable <- function (x, decreasing = FALSE, ...) {
    x$truthTable <- x$truthTable[order(x$truthTable$Consistency,decreasing=decreasing),]
    x
}
