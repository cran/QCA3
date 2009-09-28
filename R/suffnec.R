suffnec <- function(x,use=c("complete","pairwise")){
  
  consistency_fn <- function(x,y,alternative=c("xley","ylex"),...){
    ## helper function
    allvalues <- !is.na(x) & !is.na(y)
    x<-x[allvalues]
    y<-y[allvalues]
    Min <- pmin(x,y)
    alternative <- match.arg(alternative)
    Sum <- switch(
                  alternative,
                  xley=sum(x),
                  ylex=sum(y),
                  )
    ans <- sum(Min)/Sum
    return(ans)
  }
  
  minmax <- range(x,na.rm=TRUE)
  if (minmax[1] < 0 || minmax[2] > 1) stop("All the values of 'x' should be range from 0 to 1.")
  use <- match.arg(use)
  if (use=="complete") x= na.exclude(x)
  Nvar <- ncol(x)
  ans <- matrix(numeric(0),nrow=Nvar,ncol=Nvar)
  index <- t(combn(Nvar,2)) # the first col is the column index
  nindex <- nrow(index)
  
  for (i in 1:nindex) {
    rindex <- index[i,][2]
    cindex <- index[i,][1] ## row change fast -> fill lower matrix first.
    ans[rindex,cindex] <- consistency_fn(x=x[,rindex],y=x[,cindex],"xley") 
    ## sufficient condiction
    ans[cindex,rindex] <- consistency_fn(x=x[,rindex],y=x[,cindex],"ylex") 
    ## necessary condiction
  }
  dimnames(ans)= list(X=names(x),Y=names(x))
  diag(ans) <- 1
  ans2 <- t(ans)
  dimnames(ans2)= list(X=names(x),Y=names(x))
  result <- list(suff=ans,nec=ans2)
  class(result) <- "suffnec"
  result
}

print.suffnec <- function(x,digits=3,...)
{
  x<-unclass(x)
  cat("\nNecessity Scores Matrix:\n'X is necessary condition of Y'\n")
  print(x$nec,digits=digits,na.print=" ",quote = FALSE,...)
  cat("\nSufficiency Scores Matrix:\n'X is sufficient condition of Y'\n")
  print(x$suff,digits=digits,na.print=" ",quote = FALSE,...)
}

#suffnec(nssf[,c("q43d","q43e","q43f","q43g")])
