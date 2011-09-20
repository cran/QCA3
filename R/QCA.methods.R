## this is part of the QCA3 project
## by Ronggui Huang 2009-2010

excludeCSA <- function(object,csa){
  call <- match.call()
  nlevels <- object$nlevels
  conditions <- names(object$explained)
  superSets1 <- apply(object$explained,1,QCA3:::superSet,nlevels=nlevels)
  dim(superSets1) <- NULL
  superSets1 <- unique(superSets1)
  superSets0 <-  apply(QCA3:::id2Implicant(object$idExclude,nlevels),
                       1,QCA3:::superSet,nlevels=nlevels)
  dim(superSets0) <- NULL
  superSets0 <- unique(superSets0)
  superSetsCSA <- apply(csa$solutions[[1]],1, QCA3:::superSet,nlevels=nlevels)
  dim(superSetsCSA) <- NULL
  superSetsCSA <- unique(superSetsCSA)
  primesId <- sort(setdiff(superSets1, unique(c(superSets0,superSetsCSA))))
  primesId <- QCA3:::ereduce1(primesId, nlevels = nlevels)
  primeImplicants <- QCA3:::id2Implicant(primesId, nlevels = nlevels, names = conditions)
  PIChart <- QCA3:::PIChart(primeImplicants, object$explained)
  sl <- QCA3:::solvePIChart(PIChart)
  solutions <- apply(sl, 2, function(idx) primeImplicants[idx,])
  commonSolutions <- apply(sl, 1, function(idx) {
    if (length(id <- unique(idx)) == 1)
      id
  })
  ans <- list(solutions = solutions, commonSolutions = commonSolutions,
              solutionsIDX = sl, primeImplicants = primeImplicants,
              truthTable = object$truthTable, explained = object$explained,
              idExclude = object$idExclude,
              nlevels = nlevels, PIChart = PIChart,call=call)
  class(ans) <- c("QCA","noCSA")
  ans
}

primeImplicants <- function(object,traditional=TRUE){
    ## extract the prime implicants and print it in a pretty way
    nlevels <- object$nlevels
    var_names <- names(object$primeImplicants)
    PIs <- apply(object$primeImplicants, 1, QCA3:::toString, traditional = traditional,
                 nlevels = nlevels, name = var_names)
    PI <- paste(PIs, collapse = " + ")
    writeLines(strwrap(PI))
}
