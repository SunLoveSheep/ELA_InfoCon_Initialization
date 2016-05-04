#This .R file contains some functions to convert data frame, train some single models, and crossvalidation.
library("nnet")
library("neuralnet")
library("clusterSim")
library("caret")
library("kernlab")
library("RSNNS")
library("rpart")

set.seed(13)
ColumnNames <- c('hmax','epsilon_s','m0','epsilon_05')

#convert to numerical variables
Numericalization <- function(input, col_names){
  output <- data.frame(input)
  for (i in col_names) {
    output[,i] <- as.double(as.character(input[,i]))
  }
  
  return(output)
}

#normalization based on trainset's max, min:
Normalization_trainsetBased <- function(train_data, test_data, col_names, test_data_nrows = nrow(test_data))
{
  #convert to numerical format
  train_data <- Numericalization(train_data, col_names)
  test_data <- Numericalization(test_data, col_names)
  
  for (i in col_names) {
    MaxValue <- max(as.double(as.character(train_data[,i])))
    MinValue <- min(as.double(as.character(train_data[,i])))
    
    for (j in 1:test_data_nrows ) {
      test_data[j,i] <- (test_data[j,i]-MinValue)/(MaxValue-MinValue)
    }
  }
  
  return(test_data)
}

#train_data and test_data should be in numerical form
TrainNNET <- function(train_data, train_category, test_data, test_category)
{
  #Train network
  trainset_nnet <- cbind(train_data, train_category)
  colnames(trainset_nnet) <- c('hmax','epsilon_s','m0','epsilon_05','ftype')
  nnet_ftype <- nnet(ftype ~ hmax + epsilon_s + m0 + epsilon_05, trainset_nnet, size = 20, maxit = 10000)
  #View(table(predict(nnet_ftype,trainset_nnet,type="class"),trainset_nnet$ftype))
  
  #Test
  testresult_nnet <- data.frame(predict(nnet_ftype, test_data, type = "class"))
  predict_compare <- cbind(testresult_nnet, test_category)
  colnames(predict_compare) <- c('predict results','targets')
  
  #return(predict_compare)
  #return accuracy:
  accuracy <- (as.double(sum(predict_compare[,1] == predict_compare[,2]))/as.double(nrow(predict_compare)))
  return(accuracy)
}

#training using neuralnet, all input should be in numerical form (double)
TrainNeuralNet <- function(train_data, train_category, test_data, test_category)
{
  set.seed(13)
  trainset_neuralnet <- cbind(train_data,train_category)
  colnames(trainset_neuralnet) <- c('hmax','epsilon_s','m0','epsilon_05','cec14ftype')
  neuralnet_ftype <- neuralnet(cec14ftype ~ hmax + epsilon_s + m0 + epsilon_05, data = trainset_neuralnet, hidden = c(10,4), threshold = 0.001, linear.output = TRUE)
  
  testresult_neuralnet <- cbind(compute(neuralnet_ftype, test_data)$net.result, test_category)
  colnames(testresult_neuralnet) <- c('predict results','targets')
  View(testresult_neuralnet)
}

#Random Forest learning:
library(randomForest)
RF_function <- function(traindata, traintarget)
{
  rf_train = randomForest(x = traindata, y = traintarget, importance = TRUE, ntree = 101, proximity = TRUE)
  rf_confuction <- rf_train$confusion
  
  return(rf_confuction)
}

#k means clustering:
Trainkmeans <- function(input_data)
{
  kmcl <- kmeans(input_data, 2)
  plot(input_data, col = kmcl$cluster)
  
  return(kmcl)
}

#----------------------------------------------
#cross validation:
train.control <- trainControl( 
  method = "repeatedcv" 
  , number = 5 
  , repeats = 5 
  , verboseIter = T 
  , returnData = T 
  , savePredictions = T 
) 

tune.grid <- expand.grid( 
  #-----------------------
  #parameters for nnet:
  #size = c(5,10,15,20), #number of neurons
  #decay = 2^(-3:1) #decay rate
  
  #-----------------------
  #parameters for svmRBF:
  #sigma = c(0.2,0.5,0.8)
  #,C = c(0.5,1,1.5,2)
  
  #-----------------------
  #parameters for rsnns multi-layer perceptron:
  #layer1 = c(10,15,20),
  #layer2 = c(2,4,8),
  #layer3 = 0,
  #decay = 2^(-3:0)
  
  #-----------------------
  #parameters for RBF Network
  #size = c(3,4,5,6)
  
  #-----------------------
  #parameters for rpart CART
  cp = c(0.01,0.05,0.1)
) 

nnet.train <- train( 
  x = CombinedTrainset[,1:4] 
  , y = CombinedTrainset[,5] 
  , method = "rpart" #nnet(80%), svmRadial(83%+), mlpWeightDecayML(60%), rbfs(80%), rpart(CART)
  , preProcess = c("center","scale")  
  , metric = "Accuracy" 
  , trControl = train.control 
  , tuneGrid = tune.grid 
) 
nnet.train 
plot(nnet.train)
#-------------------------------------------------
