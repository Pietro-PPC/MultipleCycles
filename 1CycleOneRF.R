#Pedigree Breeding Method using GEBV Selections

#load packages
library(AlphaSimR)
library(rrBLUP)
library(caret)
library(ranger)
library(tidyverse)
library(e1071)
library(randomForest)

## Declare functions

source("declare_functions.R")

#Variables used inside functions
popList <- listenv::listenv()
EBV <- listenv::listenv()

corMat <- listenv::listenv()
allelesMat <- listenv::listenv()

#Create Results Matrices
gvMatC1 <- matrix(nrow=10, ncol=1)
corMat$C1 <- matrix(nrow=7, ncol=1)
varMatC1 <- matrix(nrow=9, ncol=1)


#establish simulation parameters
genMap <- readRDS("genMapSNPs.RData")
haplotypes <- readRDS("haplotypesSNPs.RData")

founderPop = newMapPop(genMap, 
                       haplotypes, 
                       inbred = FALSE, 
                       ploidy = 2L)

SP <- SimParam$new(founderPop)
SP$addTraitAEG(10, mean=8.8)
SP$setVarE(h2=0.25)

#INITIAL TRAINING POP

## randomly cross 200 parents 
Parents = newPop(founderPop)
TopParents = selectInd(Parents, 10, top=TRUE)

popList$F1 = randCross(TopParents, 200, nProgeny=3)

## self and bulk popList$F1 to form popList$F2 ##

popList$F2 = self(popList$F1)
popList$F2 = setPheno(popList$F2)

source("RF_F2data.R")
print("ran RF_F2data.R C1_1")

##set EBV using RF model##
setEbvRf("C1", popList, EBV, "F2", corMat, 1)

popList$NP = selectInd(popList$F2, 10, use="ebv", top=TRUE)

varMatC1[2,] = varG(popList$NP)
gvMatC1[2,] <- mean(gv(popList$NP))

setAllelesMat(allelesMat, "NP", popList)

#start new cycle

##start with 200 random crosses

popList$F1 = randCross(popList$NP, 200) 

varMatC1[3,] = varG(popList$F1)
gvMatC1[3,] <- mean(gv(popList$F1))

setAllelesMat(allelesMat, "F1", popList)

## self and bulk popList$F1 to form popList$F2 ##

popList$F2 = self(popList$F1, nProgeny=10)
varMatC1[4,] = varG(popList$F2)
gvMatC1[4,] <- mean(gv(popList$F2))

source("RF_F2data.R")

setAllelesMat(allelesMat, "F2", popList)

##set EBV using RF model##
setEbvRf("C1", popList, EBV, "F2", corMat, 2)

SelectParents = source("SelectParentsF2.R")

## select top individuals from popList$F2 bulk  to form popList$F3 ##

TopFamF2 = selectFam(popList$F2, 10, use="pheno", top=TRUE) 
SelectionsF2 = selectWithinFam(TopFamF2, 5, use="ebv", top=TRUE)

popList$F3 = self(SelectionsF2)
popList$F3 = setPheno(popList$F3)

varMatC1[5,] = varG(popList$F3)
gvMatC1[5,] <- mean(gv(popList$F3))

setAllelesMat(allelesMat, "F3", popList)

##set EBV using RF model##
setEbvRf("C1", popList, EBV, "F3", corMat, 3)

##select top within familiy from popList$F3 to form popList$F4 ##

TopFamF3 = selectFam(popList$F3,5,use="pheno", top=TRUE) 
SelectionsF3 = selectWithinFam(TopFamF3, 3, use="ebv", top=TRUE)

popList$F4 = self(SelectionsF3)
popList$F4 = setPheno(popList$F4)

varMatC1[6,] = varG(popList$F4)
gvMatC1[6,] <- mean(gv(popList$F4))

setAllelesMat(allelesMat, "F4", popList)

##set EBV using RF model##
setEbvRf("C1", popList, EBV, "F4", corMat, 4)

## select top families from popList$F4 to form popList$F5 ##

SelectionsF4 = selectFam(popList$F4, 4, use="ebv", top=TRUE)
popList$F5 = self(SelectionsF4)

varMatC1[7,]= varG(popList$F5)
gvMatC1[7,] <- mean(gv(popList$F5))

setAllelesMat(allelesMat, "F5", popList)

##set EBV using RF model##
setEbvRf("C1", popList, EBV, "F5", corMat, 5)

## select top popList$F5 families for preliminary yield trial ##

SelectionsF5 = selectFam(popList$F5, 3, use="ebv", top=TRUE) 
popList$PYT = self(SelectionsF5)
varMatC1[8,] = varG(popList$PYT)
gvMatC1[8,] <- mean(gv(popList$PYT))

setAllelesMat(allelesMat, "PYT", popList)

##set EBV using RF model##
setEbvRf("C1", popList, EBV, "PYT", corMat, 6)

## select top families from popList$PYT for popList$AYT ##

SelectionsPYT = selectFam(popList$PYT,  1, use="ebv", reps=5, top=TRUE) 
popList$AYT = self(SelectionsPYT)
varMatC1[9,] = varG(popList$AYT)
gvMatC1[9,] <- mean(gv(popList$AYT))

setAllelesMat(allelesMat, "AYT", popList)

##set EBV using RF model##
setEbvRf("C1", popList, EBV, "AYT", corMat, 7)

## select top plants to form variety ##
VarietySel = selectInd(popList$AYT, 1, use="ebv", top=TRUE)
Variety = self(VarietySel)
gvMatC1[10,] <- mean(gv(Variety))

setAllelesMat(allelesMat, "Variety", popList)

allelesMat$C1 <- rbind(allelesMat$NP, allelesMat$F1, allelesMat$F2, allelesMat$F3, allelesMat$F4, allelesMat$F5, allelesMat$PYT, allelesMat$AYT, allelesMat$Variety)

###collect bvs and ebvs###

bvebv <- listenv::listenv()

getBvebv(bvebv, "a", popList, "NP")
getBvebv(bvebv, "b", popList, "F2")
getBvebv(bvebv, "c", popList, "F3")
getBvebv(bvebv, "d", popList, "F4")
getBvebv(bvebv, "e", popList, "F5")
getBvebv(bvebv, "f", popList, "PYT")
getBvebv(bvebv, "g", popList, "AYT")

bv_ebvC1 <- rbind(bvebv$a, bvebv$b, bvebv$c, bvebv$d, bvebv$e, bvebv$f, bvebv$g)

#write files - naming convention: "model_trainingSet_descriptor_populationType_trait.csv"

source("1CycleTwo_RF.R")

#write files - naming convention: "model_trainingSet_descriptor_populationType_trait.csv"
