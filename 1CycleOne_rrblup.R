#Pedigree Breeding Method using GEBV Selections

#load packages
library(AlphaSimR)
library(rrBLUP)

#Create Results Matrices

gvMatC1 <- matrix(nrow=10, ncol=1)
corMatC1 <- matrix(nrow=7, ncol=1)
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

F1 = randCross(TopParents, 200, nProgeny=3)

## self and bulk F1 to form F2 ##

F2 = self(F1)
F2 = setPheno(F2)

source("rrblup_F2data.R")

M_F2 <-pullSegSiteGeno(F2)
G_F2 = M_F2-1
EBVF2 <- G_F2 %*% markerEffects

F2@ebv <- as.matrix(EBVF2)
corMatC1[1,] = cor(bv(F2), ebv(F2))

newParents = selectInd(F2, 10, use="ebv", top=TRUE)
varMatC1[2,] = varG(newParents)
gvMatC1[2,] <- mean(gv(newParents))

allelesMatNP <- pullSegSiteHaplo(newParents)
Gen <- as.data.frame(rep("NP", times=nInd(newParents)))
colnames(Gen) <- "Gen"
allelesMatNP <- cbind(Gen, allelesMatNP)

#start new cycle

##start with 200 random crosses

F1 = randCross(newParents, 200) 
varMatC1[3,] = varG(F1)
gvMatC1[3,] <- mean(gv(F1))

allelesMatF1 <- pullSegSiteHaplo(F1)
Gen <- as.data.frame(rep("F1", times=nInd(F1)))
colnames(Gen) <- "Gen"
allelesMatF1 <- cbind(Gen, allelesMatF1)

## self and bulk F1 to form F2 ##

F2 = self(F1, nProgeny=10) 
varMatC1[4,] = varG(F2)
gvMatC1[4,] <- mean(gv(F2))

source("rrblup_F2data.R")

allelesMatF2 <- pullSegSiteHaplo(F2)
Gen <- as.data.frame(rep("F2", times=nInd(F2)))
colnames(Gen) <- "Gen"
allelesMatF2 <- cbind(Gen, allelesMatF2)

##set EBV using BLUP model##
M_F2 <-pullSegSiteGeno(F2)
G_F2 = M_F2-1
EBVF2 <- G_F2 %*% markerEffects

F2@ebv <- as.matrix(EBVF2)
corMatC1[2,] = cor(bv(F2), ebv(F2))

SelectParents = source("SelectParentsF2.R")

## select top individuals from F2 bulk  to form F3 ##

TopFamF2 = selectFam(F2, 10, use="pheno", top=TRUE) 
SelectionsF2 = selectWithinFam(TopFamF2, 5, use="ebv", top=TRUE)

F3 = self(SelectionsF2)
F3 = setPheno(F3)
varMatC1[5,] = varG(F3)
gvMatC1[5,] <- mean(gv(F3))

allelesMatF3 <- pullSegSiteHaplo(F3)
Gen <- as.data.frame(rep("F3", times=nInd(F3)))
colnames(Gen) <- "Gen"
allelesMatF3 <- cbind(Gen, allelesMatF3)



##set EBV using BLUP model##
M_F3 <-pullSegSiteGeno(F3)
G_F3 = M_F3-1
EBVF3 <- G_F3 %*% markerEffects

F3@ebv <- as.matrix(EBVF3)
corMatC1[3,] = cor(bv(F3),ebv(F3))

##select top within familiy from F3 to form F4 ##

TopFamF3 = selectFam(F3,5,use="pheno", top=TRUE) 
SelectionsF3 = selectWithinFam(TopFamF3, 3, use="ebv", top=TRUE)

F4 = self(SelectionsF3)
F4 = setPheno(F4)
varMatC1[6,] = varG(F4)
gvMatC1[6,] <- mean(gv(F4))

allelesMatF4 <- pullSegSiteHaplo(F4)
Gen <- as.data.frame(rep("F4", times=nInd(F4)))
colnames(Gen) <- "Gen"
allelesMatF4 <- cbind(Gen, allelesMatF4)


##set EBV using BLUP model##
M_F4 <-pullSegSiteGeno(F4)
G_F4 = M_F4-1
EBVF4 <- G_F4 %*% markerEffects

F4@ebv <- as.matrix(EBVF4)
corMatC1[4,] = cor(bv(F4),ebv(F4))

## select top families from F4 to form F5 ##

SelectionsF4 = selectFam(F4, 4, use="ebv", top=TRUE)
F5 = self(SelectionsF4)
varMatC1[7,]= varG(F5)
gvMatC1[7,] <- mean(gv(F5))

allelesMatF5 <- pullSegSiteHaplo(F5)
Gen <- as.data.frame(rep("F5", times=nInd(F5)))
colnames(Gen) <- "Gen"
allelesMatF5 <- cbind(Gen, allelesMatF5)

#continue pipeline

##set EBV using BLUP model##
M_F5 <-pullSegSiteGeno(F5)
G_F5 = M_F5-1
EBVF5 <- G_F5 %*% markerEffects

F5@ebv <- as.matrix(EBVF5)
corMatC1[5,] = cor(bv(F5),ebv(F5))

## select top F5 families for preliminary yield trial ##

SelectionsF5 = selectFam(F5, 3, use="ebv", top=TRUE) 
PYT = self(SelectionsF5)
varMatC1[8,] = varG(PYT)
gvMatC1[8,] <- mean(gv(PYT))

allelesMatPYT <- pullSegSiteHaplo(PYT)
Gen <- as.data.frame(rep("PYT", times=nInd(PYT)))
colnames(Gen) <- "Gen"
allelesMatPYT <- cbind(Gen, allelesMatPYT)



##set EBV using BLUP model##
M_PYT <-pullSegSiteGeno(PYT)
G_PYT = M_PYT-1
EBVPYT <- G_PYT %*% markerEffects

PYT@ebv <- as.matrix(EBVPYT)
corMatC1[6,] = cor(bv(PYT),ebv(PYT))

## select top families from PYT for AYT ##

SelectionsPYT = selectFam(PYT,  1, use="ebv", reps=5, top=TRUE) 
AYT = self(SelectionsPYT)
varMatC1[9,] = varG(AYT)
gvMatC1[9,] <- mean(gv(AYT))

allelesMatAYT <- pullSegSiteHaplo(AYT)
Gen <- as.data.frame(rep("AYT", times=nInd(AYT)))
colnames(Gen) <- "Gen"
allelesMatAYT <- cbind(Gen, allelesMatAYT)



##set EBV using BLUP model##
M_AYT <-pullSegSiteGeno(AYT)
G_AYT = M_AYT-1
EBVAYT <- G_AYT %*% markerEffects

AYT@ebv <- as.matrix(EBVAYT)
corMatC1[7,] = cor(bv(AYT),ebv(AYT))

## select top plants to form variety ##
VarietySel = selectInd(AYT, 1, use="ebv", top=TRUE)
Variety = self(VarietySel)
gvMatC1[10,] <- mean(gv(Variety))

allelesMatVar <- pullSegSiteHaplo(Variety)
Gen <- as.data.frame(rep("Variety", times=nInd(Variety)))
colnames(Gen) <- "Gen"
allelesMatVar <- cbind(Gen, allelesMatVar)

allelesMatC1 <- rbind(allelesMatNP, allelesMatF1, allelesMatF2, allelesMatF3, allelesMatF4, allelesMatF5, allelesMatPYT, allelesMatAYT, allelesMatVar)


###collect bvs and ebvs###

bvebv <- cbind(bv(newParents), ebv(newParents))
Gen <- as.data.frame(rep("NP", times=nInd(newParents)))
bvebv <- cbind(Gen, bvebv)
colnames(bvebv) <- c("Gen","bv","ebv")

bvebv1 <- cbind(bv(F2), ebv(F2))
Gen <- as.data.frame(rep("F2", times=nInd(F2)))
bvebv1 <- cbind(Gen, bvebv1)
colnames(bvebv1) <- c("Gen","bv","ebv")

bvebv2 <- cbind(bv(F3), ebv(F3))
Gen <- as.data.frame(rep("F3", times=nInd(F3)))
bvebv2 <- cbind(Gen, bvebv2)
colnames(bvebv2) <- c("Gen","bv","ebv")

bvebv3 <- cbind(bv(F4), ebv(F4))
Gen <- as.data.frame(rep("F4", times=nInd(F4)))
bvebv3 <- cbind(Gen, bvebv3)
colnames(bvebv3) <- c("Gen","bv","ebv")

bvebv4 <- cbind(bv(F5), ebv(F5))
Gen <- as.data.frame(rep("F5", times=nInd(F5)))
bvebv4 <- cbind(Gen, bvebv4)
colnames(bvebv4) <- c("Gen","bv","ebv")

bvebv5 <- cbind(bv(PYT), ebv(PYT))
Gen <- as.data.frame(rep("PYT", times=nInd(PYT)))
bvebv5 <- cbind(Gen, bvebv5)
colnames(bvebv5) <- c("Gen","bv","ebv")

bvebv6 <- cbind(bv(AYT), ebv(AYT))
Gen <- as.data.frame(rep("AYT", times=nInd(AYT)))
bvebv6 <- cbind(Gen, bvebv6)
colnames(bvebv6) <- c("Gen","bv","ebv")


bv_ebvC1 <- rbind(bvebv,bvebv1,bvebv2,bvebv3,bvebv4,bvebv5,bvebv6)


#write files - naming convention: "model_trainingSet_descriptor_populationType_trait.csv"



source("1CycleThree_rrblup.R")

