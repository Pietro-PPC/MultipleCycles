
TrainingGeno <- pullSegSiteGeno(F2)
TrainingPheno <- pheno(F2)

pheno <- as.data.frame(TrainingPheno)
geno <- as.data.frame(TrainingGeno)
trainingdata <- cbind(pheno, geno)
colnames(trainingdata) <- paste("ID",1:ncol(trainingdata), sep="") 


## create TRN TST split ###
train_index <- sample(1:nrow(trainingdata), 0.75 * nrow(trainingdata))
trainingset <- trainingdata[train_index, ]
testingset <- trainingdata[-train_index, ]


## create cross validation strategy ##
control <- trainControl(method='repeatedcv', 
                        number=10, ##will test 10 different values for mtry (number of variables for splitting) ##
                        repeats=3,
                        search = "random")

##build model##
cl <- makePSOCKcluster(15)
print('started cluster')
registerDoParallel(cl)
print('registered cluster')

rf_fit = train(ID1 ~ .,
               data = trainingset,
               method = "parRF",
               tuneLength = 10,
               trControl=control) ## search a random tuning grid ##

print('Ran in parallel')
stopCluster(cl)
print('Closed cluster')

### This command takes about 90 minutes in an compute canada interactive session ###
