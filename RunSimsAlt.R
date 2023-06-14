## define variables ##

nModels = 7
nReps = 1
nGen = 10
nVar = 9


geneticvalues <- listenv::listenv()
correlations <- listenv::listenv()
variances <- listenv::listenv()
alleles <- listenv::listenv()
bvebv <- listenv::listenv()

## establish empty matrices to hold outputs for Selfing and Recombination Population ##
for (cycle in c("C1","C2","C3")){
  geneticvalues[[cycle]] <- matrix(nrow=nGen, ncol=nReps)  
  correlations[[cycle]] <- matrix(nrow=nModels, ncol=nReps)
  variances[[cycle]] <- matrix(nrow=nVar,ncol=nReps)
  alleles[[cycle]] <- vector("list", length = nReps)
  bvebv[[cycle]] <- vector("list", length = nReps)
}

## Run repeat loop to run reps ##

for(i in 1:nReps){
  source("CyclesRF.R") ##Source the SCript for the SCenario you would like to run##
  
  for (cycle in c("C1","C2","C3")){
    geneticvalues[[cycle]][,i] <- gvMat[[cycle]]  
    correlations[[cycle]][,i] <- corMat[[cycle]]
    variances[[cycle]][,i] <- varMat[[cycle]]
    alleles[[cycle]][[i]] <- allelesMat[[cycle]]
    bvebv[[cycle]][[i]] <- bv_ebv[[cycle]]
  }
}
  
##create data frames and label##
for (cycle in c("C1", "C2", "C3")){
  geneticvalues[[cycle]] <- as.data.frame(geneticvalues[[cycle]])
  colnames(geneticvalues[[cycle]]) <- 1:nReps

  if (cycle == "C1")
    gain <- as.data.frame(geneticvalues[[cycle]][10,] - geneticvalues[[cycle]][2,])
  else 
    gain <- as.data.frame(geneticvalues[[cycle]][10,] - geneticvalues[[cycle]][3,])
  colnames(gain) <- 1:nReps

  Allgeneticvalues[[cycle]] <- as.data.frame(rbind(geneticvalues[[cycle]], gain))
  rownames(Allgeneticvalues[[cycle]]) <- c("PrevCycPYT","NewParents","F1","F2","F3","F4","F5","PYT","AYT","Variety","meanGV")
  colnames(Allgeneticvalues[[cycle]]) <- c(1:nReps)

  correlations[[cycle]] <- as.data.frame(correlations[[cycle]])
  rownames(correlations[[cycle]]) <- c("NewParents","F2","F3","F4","F5","PYT","AYT")
  colnames(correlations[[cycle]]) <- c(1:nReps)

  variances[[cycle]] <- as.data.frame(variances[[cycle]])
  colnames(variances[[cycle]]) <- c(1:nReps)
  rownames(variances[[cycle]]) <- c("PrevCycPYT", "newParents","F1","F2", "F3","F4", "F5", "PYT","AYT")

  ##write files

  write.csv(Allgeneticvalues[[cycle]], paste("1", cycle, "_rrblup_rd_gvs_snp_yield.csv", sep=""))
  write.csv(correlations[[cycle]], paste("1", cycle, "_rrblup_rd_cors_snp_yield.csv", sep=""))
  write.csv(variances[[cycle]], paste("1", cycle, "_rrblup_rd_vars_snp_yield.csv", sep=""))
  saveRDS(alleles[[cycle]], file=paste("1", cycle, "rrblup_rd_alleles_snp_yield.rds", sep=""))
  saveRDS(bvebv[[cycle]], file=paste("1", cycle, "rrblup_rd_bvebv_snp_yield.rds", sep=""))
}
