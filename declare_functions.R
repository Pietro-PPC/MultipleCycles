## Function declarations

setEbvRf <- function(cycle, popList, EBV, pop_name, corMat, corMat_ind){
    # Sets EBV for cycles that use Random Forest
    ## popList is a list of populations be F2, F3, ...
    ## EBV is a list of EBVs
    ## pop_name can be "F2", "F3", ...
    ## corMat is a list of matrices "C1", "C2", ...
    ## corMat_ind can be 2, 3, ...
    M = as.data.frame(pullSegSiteGeno(popList[[pop_name]]))
    colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
    EBV[[pop_name]] <- as.numeric(predict(rf_fit, M)) # rf_fit is the result from random forest

    popList[[pop_name]]@ebv <- as.matrix(EBV[[pop_name]])
    corMat[[cycle]][corMat_ind,] = cor(bv(popList[[pop_name]]), ebv(popList[[pop_name]]))
}

# Gets bvebv and assigns to list
getBvebv <- function(bvebv, key, popList, pop_name){
    bvebv[[key]] <- cbind(bv(popList[[pop_name]]), ebv(popList[[pop_name]]))
    Gen <- as.data.frame(rep(pop_name, times=nInd(popList[[pop_name]])))
    bvebv[[key]] <- cbind(Gen, bvebv[[key]])
    colnames(bvebv[[key]]) <- c("Gen","bv","ebv")
}

# Adds alleles matrix data for population pop_name
setAllelesMat <- function(mat, pop_name, popList){
    allelesMat[[pop_name]] <- pullSegSiteHaplo(popList[[pop_name]])
    Gen <- as.data.frame(rep(pop_name, times=nInd(popList[[pop_name]])))
    colnames(Gen) <- "Gen"
    allelesMat[[pop_name]] <- cbind(Gen, allelesMat[[pop_name]])
}

