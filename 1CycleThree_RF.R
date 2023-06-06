corMatC3 <- listenv::listenv()

gvMatC3 <- matrix(nrow=10, ncol=1)
corMatC3$m <- matrix(nrow=7, ncol=1)
varMatC3 <- matrix(nrow=9, ncol=1)

corMatC3$m[1] = cor(bv(newCycleSelections), ebv(newCycleSelections))

popList[["F1"]] = makeCross(newCycleSelections, crossPlan = cross, nProgeny = 5)
varMatC3[3,] = varG(popList[["F1"]])
gvMatC3[3,] <- mean(gv(popList[["F1"]]))

allelesMatF1 <- pullSegSiteHaplo(popList[["F1"]])
Gen <- as.data.frame(rep("F1", times=nInd(popList[["F1"]])))
colnames(Gen) <- "Gen"
allelesMatF1 <- cbind(Gen, allelesMatF1)

## self and bulk popList[["F1"]] to form popList[["F2"]] ##

popList[["F2"]] = self(popList[["F1"]], nProgeny = 10) 
varMatC3[4,] = varG(popList[["F2"]])
gvMatC3[4,] <- mean(gv(popList[["F2"]]))

allelesMatF2 <- pullSegSiteHaplo(popList[["F2"]])
Gen <- as.data.frame(rep("F2", times=nInd(popList[["F2"]])))
colnames(Gen) <- "Gen"
allelesMatF2 <- cbind(Gen, allelesMatF2)

source("RF_F2data.R")
print("ran RF_F2data.R C3")

M = as.data.frame(pullSegSiteGeno(popList[["F2"]]))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
EBVF2 <- as.numeric(predict(rf_fit, M))


popList[["F2"]]@ebv <- as.matrix(EBVF2)
corMatC3$m[2] = cor(bv(popList[["F2"]]), ebv(popList[["F2"]]))


rm(newCycleSelections)
source("SelectParentsF2.R")

## select top individuals from popList[["F2"]] bulk  to form popList[["F3"]] ##

TopFamF2 = selectFam(popList[["F2"]], 10, use="ebv", top=TRUE) 
SelectionsF2 = selectWithinFam(TopFamF2, 5, use="ebv", top=TRUE)

popList[["F3"]] = self(SelectionsF2)
popList[["F3"]] = setPheno(popList[["F3"]])
varMatC3[5,] = varG(popList[["F3"]])
gvMatC3[5,] <- mean(gv(popList[["F3"]]))

allelesMatF3 <- pullSegSiteHaplo(popList[["F3"]])
Gen <- as.data.frame(rep("F3", times=nInd(popList[["F3"]])))
colnames(Gen) <- "Gen"
allelesMatF3 <- cbind(Gen, allelesMatF3)

#set EBV using RF model
setEbvRf(popList, EBV, "F3", corMatC3, 3)

##select top within familiy from popList[["F3"]] to form popList[["F4"]] ##

TopFamF3 = selectFam(popList[["F3"]],5,use="pheno", top=TRUE) 
SelectionsF3 = selectWithinFam(TopFamF3, 3, use="ebv", top=TRUE)

popList[["F4"]] = self(SelectionsF3)
popList[["F4"]] = setPheno(popList[["F4"]])
varMatC3[6,] = varG(popList[["F4"]])
gvMatC3[6,] <- mean(gv(popList[["F4"]]))

allelesMatF4 <- pullSegSiteHaplo(popList[["F4"]])
Gen <- as.data.frame(rep("F4", times=nInd(popList[["F4"]])))
colnames(Gen) <- "Gen"
allelesMatF4 <- cbind(Gen, allelesMatF4)


##set EBV using BLUP model##
setEbvRf(popList, EBV, "F4", corMatC3, 4)

## select top families from popList[["F4"]] to form popList[["F5"]] ##

SelectionsF4 = selectFam(popList[["F4"]], 4, use="ebv", top=TRUE)
popList[["F5"]] = self(SelectionsF4)
varMatC3[7,]= varG(popList[["F5"]])
gvMatC3[7,] <- mean(gv(popList[["F5"]]))

allelesMatF5 <- pullSegSiteHaplo(popList[["F5"]])
Gen <- as.data.frame(rep("F5", times=nInd(popList[["F5"]])))
colnames(Gen) <- "Gen"
allelesMatF5 <- cbind(Gen, allelesMatF5)


#continue pipeline

##set EBV using BLUP model##
setEbvRf(popList, EBV, "F5", corMatC3, 5)

## select top popList[["F5"]] families for preliminary yield trial ##

SelectionsF5 = selectFam(popList[["F5"]], 3, use="ebv", top=TRUE) 
popList[["PYT"]] = self(SelectionsF5)
varMatC3[8,] = varG(popList[["PYT"]])
gvMatC3[8,] <- mean(gv(popList[["PYT"]]))

allelesMatPYT <- pullSegSiteHaplo(popList[["PYT"]])
Gen <- as.data.frame(rep("PYT", times=nInd(popList[["PYT"]])))
colnames(Gen) <- "Gen"
allelesMatPYT <- cbind(Gen, allelesMatPYT)

##set EBV using BLUP model##
setEbvRf(popList, EBV, "PYT", corMatC3, 6)

## select top families from popList[["PYT"]] for popList[["AYT"]] ##

SelectionsPYT = selectFam(popList[["PYT"]],  1, use="ebv", reps=5, top=TRUE) 
popList[["AYT"]] = self(SelectionsPYT)
varMatC3[9,] = varG(popList[["AYT"]])
gvMatC3[9,] <- mean(gv(popList[["AYT"]]))

allelesMatAYT <- pullSegSiteHaplo(popList[["AYT"]])
Gen <- as.data.frame(rep("AYT", times=nInd(popList[["AYT"]])))
colnames(Gen) <- "Gen"
allelesMatAYT <- cbind(Gen, allelesMatAYT)

##set EBV using BLUP model##
setEbvRf(popList, EBV, "AYT", corMatC3, 7)

## select top plants to form variety ##
VarietySel = selectInd(popList[["AYT"]], 1, use="ebv", top=TRUE)
Variety = self(VarietySel)
gvMatC3[10,] <- mean(gv(Variety))

allelesMatVar <- pullSegSiteHaplo(Variety)
Gen <- as.data.frame(rep("Variety", times=nInd(Variety)))
colnames(Gen) <- "Gen"
allelesMatVar <- cbind(Gen, allelesMatVar)

allelesMatC3 <- rbind(allelesMatNP, allelesMatF1, allelesMatF2, allelesMatF3, allelesMatF4, allelesMatF5, allelesMatPYT, allelesMatAYT, allelesMatVar)

###collect bvs and ebvs###

bvebv <- cbind(bv(newParents), ebv(newParents))
Gen <- as.data.frame(rep("NP", times=nInd(newParents)))
bvebv <- cbind(Gen, bvebv)
colnames(bvebv) <- c("Gen","bv","ebv")

bvebv1 <- cbind(bv(popList[["F2"]]), ebv(popList[["F2"]]))
Gen <- as.data.frame(rep("F2", times=nInd(popList[["F2"]])))
bvebv1 <- cbind(Gen, bvebv1)
colnames(bvebv1) <- c("Gen","bv","ebv")

bvebv2 <- cbind(bv(popList[["F3"]]), ebv(popList[["F3"]]))
Gen <- as.data.frame(rep("F3", times=nInd(popList[["F3"]])))
bvebv2 <- cbind(Gen, bvebv2)
colnames(bvebv2) <- c("Gen","bv","ebv")

bvebv3 <- cbind(bv(popList[["F4"]]), ebv(popList[["F4"]]))
Gen <- as.data.frame(rep("F4", times=nInd(popList[["F4"]])))
bvebv3 <- cbind(Gen, bvebv3)
colnames(bvebv3) <- c("Gen","bv","ebv")

bvebv4 <- cbind(bv(popList[["F5"]]), ebv(popList[["F5"]]))
Gen <- as.data.frame(rep("F5", times=nInd(popList[["F5"]])))
bvebv4 <- cbind(Gen, bvebv4)
colnames(bvebv4) <- c("Gen","bv","ebv")

bvebv5 <- cbind(bv(popList[["PYT"]]), ebv(popList[["PYT"]]))
Gen <- as.data.frame(rep("PYT", times=nInd(popList[["PYT"]])))
bvebv5 <- cbind(Gen, bvebv5)
colnames(bvebv5) <- c("Gen","bv","ebv")

bvebv6 <- cbind(bv(popList[["AYT"]]), ebv(popList[["AYT"]]))
Gen <- as.data.frame(rep("AYT", times=nInd(popList[["AYT"]])))
bvebv6 <- cbind(Gen, bvebv6)
colnames(bvebv6) <- c("Gen","bv","ebv")


bv_ebvC2 <- rbind(bvebv,bvebv1,bvebv2,bvebv3,bvebv4,bvebv5,bvebv6)



#write files - naming convention: "model_trainingSet_descriptor_populationType_trait.csv"



