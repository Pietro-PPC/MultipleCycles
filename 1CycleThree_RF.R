gvMatC3 <- matrix(nrow=10, ncol=1)
corMat$C3 <- matrix(nrow=7, ncol=1)
varMatC3 <- matrix(nrow=9, ncol=1)

corMat$C3[1] = cor(bv(newCycleSelections), ebv(newCycleSelections))

popList$F1 = makeCross(newCycleSelections, crossPlan = cross, nProgeny = 5)
varMatC3[3,] = varG(popList$F1)
gvMatC3[3,] <- mean(gv(popList$F1))

setAllelesMat(allelesMat, "F1", popList)

## self and bulk popList$F1 to form popList$F2 ##

popList$F2 = self(popList$F1, nProgeny = 10) 
varMatC3[4,] = varG(popList$F2)
gvMatC3[4,] <- mean(gv(popList$F2))

setAllelesMat(allelesMat, "F2", popList)

source("RF_F2data.R")
print("ran RF_F2data.R C3")

M = as.data.frame(pullSegSiteGeno(popList$F2))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
EBVF2 <- as.numeric(predict(rf_fit, M))


popList$F2@ebv <- as.matrix(EBVF2)
corMat$C3[2] = cor(bv(popList$F2), ebv(popList$F2))


rm(newCycleSelections)
source("SelectParentsF2.R")

## select top individuals from popList$F2 bulk  to form popList$F3 ##

TopFamF2 = selectFam(popList$F2, 10, use="ebv", top=TRUE) 
SelectionsF2 = selectWithinFam(TopFamF2, 5, use="ebv", top=TRUE)

popList$F3 = self(SelectionsF2)
popList$F3 = setPheno(popList$F3)
varMatC3[5,] = varG(popList$F3)
gvMatC3[5,] <- mean(gv(popList$F3))

setAllelesMat(allelesMat, "F3", popList)

#set EBV using RF model
setEbvRf("C3", popList, EBV, "F3", corMat, 3)

##select top within familiy from popList$F3 to form popList$F4 ##

TopFamF3 = selectFam(popList$F3,5,use="pheno", top=TRUE) 
SelectionsF3 = selectWithinFam(TopFamF3, 3, use="ebv", top=TRUE)

popList$F4 = self(SelectionsF3)
popList$F4 = setPheno(popList$F4)
varMatC3[6,] = varG(popList$F4)
gvMatC3[6,] <- mean(gv(popList$F4))

setAllelesMat(allelesMat, "F4", popList)

##set EBV using BLUP model##
setEbvRf("C3", popList, EBV, "F4", corMat, 4)

## select top families from popList$F4 to form popList$F5 ##

SelectionsF4 = selectFam(popList$F4, 4, use="ebv", top=TRUE)
popList$F5 = self(SelectionsF4)
varMatC3[7,]= varG(popList$F5)
gvMatC3[7,] <- mean(gv(popList$F5))

setAllelesMat(allelesMat, "F5", popList)

#continue pipeline

##set EBV using BLUP model##
setEbvRf("C3", popList, EBV, "F5", corMat, 5)

## select top popList$F5 families for preliminary yield trial ##

SelectionsF5 = selectFam(popList$F5, 3, use="ebv", top=TRUE) 
popList$PYT = self(SelectionsF5)
varMatC3[8,] = varG(popList$PYT)
gvMatC3[8,] <- mean(gv(popList$PYT))

setAllelesMat(allelesMat, "PYT", popList)

##set EBV using BLUP model##
setEbvRf("C3", popList, EBV, "PYT", corMat, 6)

## select top families from popList$PYT for popList$AYT ##

SelectionsPYT = selectFam(popList$PYT,  1, use="ebv", reps=5, top=TRUE) 
popList$AYT = self(SelectionsPYT)
varMatC3[9,] = varG(popList$AYT)
gvMatC3[9,] <- mean(gv(popList$AYT))

setAllelesMat(allelesMat, "AYT", popList)

##set EBV using BLUP model##
setEbvRf("C3", popList, EBV, "AYT", corMat, 7)

## select top plants to form variety ##
VarietySel = selectInd(popList$AYT, 1, use="ebv", top=TRUE)
Variety = self(VarietySel)
gvMatC3[10,] <- mean(gv(Variety))

setAllelesMat(allelesMat, "Variety", popList)

allelesMat$C3 <- rbind(allelesMat$NP, allelesMat$F1, allelesMat$F2, allelesMat$F3, allelesMat$F4, allelesMat$F5, allelesMat$PYT, allelesMat$AYT, allelesMat$Var)

###collect bvs and ebvs###

getBvebv(bvebv, "a", popList, "NP")
getBvebv(bvebv, "b", popList, "F2")
getBvebv(bvebv, "c", popList, "F3")
getBvebv(bvebv, "d", popList, "F4")
getBvebv(bvebv, "e", popList, "F5")
getBvebv(bvebv, "f", popList, "PYT")
getBvebv(bvebv, "g", popList, "AYT")

bv_ebvC3 <- rbind(bvebv$a, bvebv$b, bvebv$c, bvebv$d, bvebv$e, bvebv$f, bvebv$g)

#write files - naming convention: "model_trainingSet_descriptor_populationType_trait.csv"
