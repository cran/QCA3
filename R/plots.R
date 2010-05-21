plot.QCA <- function(x,...){
    explain <- x$call$explain
    truthTable <- x$truthTable
    if (pmatch(explain,"positive",0)==1) Case <- truthTable[rownames(truthTable) %in% rownames(x$explained),]
    if (pmatch(explain,"negative",0)==1) Case <- truthTable[rownames(truthTable) %in% rownames(x$explained),]
    conditions <- names(x$explained)
    idExplained <- apply(Case[,names(x$explained)],1,QCA3:::implicant2Id,nlevel=x$nlevels)
    ids <- rep(idExplained,Case$N)
    Coverage <- apply(x$solutions[[1]],1, function(x) {
        common <- intersect(QCA3:::esubSet(x),ids)
        ## subSet or subCombination?
        ids[ids %in% common]
    }
                      )
    names(Coverage)<- paste("IM",seq_len(length(Coverage)),sep=".")
    CovList <- lapply(Coverage,as.character)
    require(venneuler) ## use venneuler instead of Vennerable
    m <- data.frame(element=unlist(CovList),
                    set=rep(names(CovList),lapply(CovList,length)))
    v <- venneuler::venneuler(m)
    venneuler:::plot.VennDiagram(v)
    ##require(Vennerable)
    ##plot(Venn(Sets=CovList),doWeights=TRUE)
}
