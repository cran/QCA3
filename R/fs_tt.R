fs_tt <- function(mydata, outcome = "", conditions = c(""),ncases_cutoff=1,consistency_cutoff=0.8,
                          complete = FALSE,show.cases = TRUE, quiet = FALSE,cases=NULL)
  ## change the name to fs_tt, fs_tt is used to construct a truthTable of class tt which is in QCA package.
  ## fs_truthTable is reserved for ASRR package.
{
  require(QCA)
  membership_cutoff=0.5
  if (outcome==""||conditions=="") stop("You must specific outcome and conditions first.")
  fulldata <- mydata[,c(outcome,conditions)]
  if (any(fulldata<0)|| any(fulldata>1)) stop("Fuzzy set score must in [0,1].")
  
  allExpress <- eval(parse(text=(sprintf("expand.grid(%s)",paste(conditions,"=1:0",sep="",collapse=",")))))
  conditions <- mydata[,conditions]
  
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
  
  score_mat <- apply(allExpress,1,function(x) getScore(x,data=conditions))
  
  allExpress$N_Cases<- apply(score_mat,2,function(x) sum(x>membership_cutoff))
  allExpress$Consistency <- apply(score_mat,2,function(x,outcome) {sum(pmin(x,outcome))/sum(x)},outcome=mydata[,outcome])
  
  allExpress$OUT <- "?"
  allExpress$OUT[allExpress$N_Cases >= ncases_cutoff & allExpress$Consistency > consistency_cutoff]<-"1"
  allExpress$OUT[allExpress$N_Cases >= ncases_cutoff & allExpress$Consistency <= consistency_cutoff]<-"0"
  allExpress$freq0 <- allExpress$freq1 <- "-"
  allExpress$freq0[allExpress$OUT=="0"] <- allExpress$N_Cases[allExpress$OUT=="0"]
  allExpress$freq1[allExpress$OUT=="1"] <- allExpress$N_Cases[allExpress$OUT=="1"]
  allExpress <- allExpress[,c(seq_len(length(conditions)),(length(conditions)+3):(length(conditions)+5),(length(conditions)+1):(length(conditions)+2))]
  ## reorder alExpress

  if (show.cases){
    if (is.null(cases)) cases <- rownames(mydata) else cases <- mydata[,cases]
    cases <- gsub(",","_",cases)
    allExpress$cases <- apply(score_mat,2,function(x) paste(cases[which( x > membership_cutoff)],sep="",collapse=","))
    if (!complete) allExpress <- allExpress[allExpress$OUT != "?",,drop=FALSE]
    tt <- list(tt=allExpress,indexes=which(allExpress$OUT!="?"),noflevels=rep(2,length(conditions)),casenames=allExpress$cases)## should use full range of name now, the bug has been fixed in eqmcc
  } else {  
    if (!complete) allExpress <- allExpress[allExpress$OUT != "?",,drop=FALSE]
    tt <- list(tt=allExpress ,indexes=which(allExpress$OUT!="?"),noflevels=rep(2,length(conditions)))
  }
  class(tt) <-"tt"
  tt
}
