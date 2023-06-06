## EBV must be a list

list_env <- listenv::listenv()

setEbvRf <- function(popList, EBV, pop_name, corMat_ind){
    ## popList is a list of populations be F2, F3, ...
    ## EBV is a list of EBVs
    ## pop_name can be "F2", "F3", ...
    ## corMat_ind can be 2, 3, ...
    M = as.data.frame(pullSegSiteGeno(popList[[pop_name]]))
    colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
    EBV[[pop_name]] <- as.numeric(predict(rf_fit, M))

    popList[[pop_name]]@ebv <- as.matrix(EBV[[pop_name]])
    corMatC1[corMat_ind,] = cor(bv(popList[[pop_name]]), ebv(popList[[pop_name]]))
}

setEbvRf(pop, EBV, "F2", 2)
setEbvRf(pop, EBV, "F3", 3)
setEbvRf(pop, EBV, "F4", 4)
setEbvRf(pop, EBV, "F5", 5)
setEbvRf(pop, EBV, "PYT", 6)
setEbvRf(pop, EBV, "AYT", 7)

M = as.data.frame(pullSegSiteGeno(F2))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
EBVF2 <- as.numeric(predict(rf_fit, M))

F2@ebv <- as.matrix(EBVF2)
corMatC1[2,] = cor(bv(F2), ebv(F2))
##
M = as.data.frame(pullSegSiteGeno(F3))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
EBVF3 <- as.numeric(predict(rf_fit, M))

F3@ebv <- as.matrix(EBVF3)
corMatC1[3,] = cor(bv(F3),ebv(F3))
##
M = as.data.frame(pullSegSiteGeno(F4))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
EBVF4 <- as.numeric(predict(rf_fit, M))

F4@ebv <- as.matrix(EBVF4)
corMatC1[4,] = cor(bv(F4),ebv(F4))
##

M = as.data.frame(pullSegSiteGeno(F5))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
EBVF5 <- as.numeric(predict(rf_fit, M))

F5@ebv <- as.matrix(EBVF5)
corMatC1[5,] = cor(bv(F5),ebv(F5))

##
M = as.data.frame(pullSegSiteGeno(PYT))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
EBVPYT <- as.numeric(predict(rf_fit, M))

PYT@ebv <- as.matrix(EBVPYT)
corMatC1[6,] = cor(bv(PYT),ebv(PYT))

##
M = as.data.frame(pullSegSiteGeno(AYT))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
EBVAYT <- as.numeric(predict(rf_fit, M))


AYT@ebv <- as.matrix(EBVAYT)
corMatC1[7,] = cor(bv(AYT),ebv(AYT))