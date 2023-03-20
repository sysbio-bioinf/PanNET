# R-version: 4.2.2

library(TunePareto) # package version: 2.5.2
library(ORION) # package version: 1.0.3
library(e1071) # package version: 1.7-13
library(GEOquery) # package version: 2.66.0       # If this is not installed, first install the package "BiocManager",
                                                  # then use BiocManager::install("GEOquery").

library(THRESHOLDerrorMinimizer) # local library, see readMe file 




###########################################################################################
# Loading the datasets GSE73338 and GSE117851 
###########################################################################################
options('download.file.method.GEOquery' = 'libcurl')

###########################################################################################
# GSE73338
###########################################################################################
GSE73338 <- getGEO(GEO = "GSE73338", GSEMatrix = TRUE, destdir=getwd())

# Reformating the list, as well as preprocessing the data 
GSE73338[["name"]] <- "GSE73338"
# For removing all samples labeled "Normal pancreas islet", as well as the single sample labeled "Functional":
sampleRM.73338 <- c(which(GSE73338[["GSE73338_series_matrix.txt.gz"]]@phenoData@data[["source_name_ch1"]] == "Normal pancreas islet"), 
                    which(GSE73338[["GSE73338_series_matrix.txt.gz"]]@phenoData@data[["source_name_ch1"]] == "Functional"))
GSE73338[["data"]] <- GSE73338[["GSE73338_series_matrix.txt.gz"]]@assayData[["exprs"]][,-sampleRM.73338]
GSE73338[["labs"]] <- GSE73338[["GSE73338_series_matrix.txt.gz"]]@phenoData@data[["source_name_ch1"]][-sampleRM.73338]
GSE73338[["samples"]] <- GSE73338[["GSE73338_series_matrix.txt.gz"]]@phenoData@data[["geo_accession"]][-sampleRM.73338]
GSE73338[["id"]] <- GSE73338[["GSE73338_series_matrix.txt.gz"]]@featureData@data[["ID"]]




###########################################################################################
# GSE117851
###########################################################################################
GSE117851 <- getGEO(GEO = "GSE117851", GSEMatrix = TRUE, destdir=getwd())

GSE117851[["name"]] <- "GSE117851"
# For remove all "genotype: ATRX" as well as "genotype: ATRX and MEN1"
sampleRM.117851 <- c(which( GSE117851[["GSE117851_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.1"]] == "genotype: ATRX"), 
                     which( GSE117851[["GSE117851_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.1"]] == "genotype: ATRX and MEN1"))
GSE117851[["data"]] <- GSE117851[["GSE117851_series_matrix.txt.gz"]]@assayData[["exprs"]][,-sampleRM.117851]
GSE117851[["labs"]] <- GSE117851[["GSE117851_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.1"]][-sampleRM.117851]
GSE117851[["samples"]] <- GSE117851[["GSE117851_series_matrix.txt.gz"]]@phenoData@data[["geo_accession"]][-sampleRM.117851]
GSE117851[["id"]] <- GSE117851[["GSE117851_series_matrix.txt.gz"]]@featureData@data[["ID"]]









###########################################################################################
# Defining the wrapper opject for calling the one-vs-one-SVM 
###########################################################################################

tunePareto.fssvm <- function(){  
  tuneParetoClassifier(   name = "fssvm",
                          classifier = function(data, labs, numFea, type = "classwise", direction = "b", scale = TRUE)
                          {
                            
                            classes <- unique(labs)
                            
                            my.labs <- rep(0,length(labs))
                            my.labs[labs == classes[1]] <- 1
                            
                            # Computing individual thresholds on each feature 
                            model <- baselearner.threshold(t(data), my.labs,type = type, direction = direction)
                            
                            # Sorting these features by their error rate 
                            score <- as.numeric(model$base.learners[,"Score"])
                            score[is.na(score)] <- 1
                            
                            # A number of numFea features with the smallest error rates are selected to use in training 
                            features <- order(score)[1:numFea]
                            
                            if(length(features)<numFea)
                            {   warning("Not enough features selected.", paste(features, collapse = " "))
                              features <- c(features, 1:numFea)[1:numFea]
                            }
                            
                            # Training of the binary SVM on the selected features
                            model <- svm(x = data[,features,drop=FALSE], y = labs, scale = scale, type = "C", kernel = "linear")
                          
                            
                            object = list(  model = model,
                                            features = features)
                            
                            return(object)
                          },
                          ,
                          classifierParamNames = c("numFea","type","direction","scale"),
                          predictor = function(object, newdata)
                          {
                            require(e1071)
                            result <- predict(object = object$model, newdata = newdata[,object$features,drop=FALSE])
                            
                            if(any(is.na(result)) | any(is.nan(result)) | any(is.infinite(result)))
                            {
                              stop("wrong prediction")
                            }
                            
                            return(result)
                          },
                          trainDataName = "data",
                          trainLabelName = "labs",
                          testDataName = "newdata",
                          modelName = "object",
                          requiredPackages = c("THRESHOLDerrorMinimizer","e1071") 
  )
}




###########################################################################################
# Function for computing the pairwise linear one-vs-one-SVMs
    # Input:    dataset - a list comprised of the following elements: 
    #                       - name: String declaring the dataset name 
    #                       - data: the dataset with the columns being the samples and 
    #                               the rows the features 
    #                       - labs: a vector of length equal to the number of samples, here 
    #                               the i-th entry is the label associated with the i-th sample
    #                       - id: a vector of length equal to the number of features, here the 
    #                             i-th entry is the name of the i-th feature
    #           thresh -  a set value in [0,1] for utilizing the ORION package computing the  
    #                     binary classification as well as the foreign classification
    # Output:   list(predMap,subc)
    #           predMap - depicts the results of the binary classification training 
    #                     using a cross validation set up, as  well as the evaluation of all samples 
    #                     using the trained classifiers
    #           subc - depicts the sensitivities of the pairwise binary classification
###########################################################################################


start.run <- function(dataset, thresh = 0.000001){
  
  #Renaming the labels of the provided dataset in order to make use of the ORION Package 
  classes <- unique(sort(dataset$labs))
  
  new.labs <- dataset$labs
  list.labs.association <- vector("list", length(classes))
  for(i in 1:length(classes))
  {
    new.labs[dataset$labs==classes[i]] <- i-1
    list.labs.association[i] <- paste(classes[i], "=", i-1)
  }
  dataset$labs <- as.numeric(new.labs)
  
  # Construction a set of partitions for a cross-validation 
  foldList <- generateCVRuns(dataset$labs, stratified = TRUE, leaveOneOut = TRUE)
  save(foldList, file = paste(dataset$name,"_foldList.RData", sep = ""))
  
  
  # Pairwise training for all binary class combinations and evaluates all samples using the trained classifiers
  predMap = predictionMap( data=t(dataset$data),
                             labels=dataset$labs,
                             foldList = foldList,
                             classifier = tunePareto.fssvm(),
                             parallel = FALSE,
                             numFea = 200)

    
  save(predMap, file = paste(dataset$name,"_predMap.RData", sep = ""))
  
  
  
  # Computes the sensitivities of the pairwise binary classification
  subc = subcascades( predictionMap = predMap,
                      sets = NULL,
                      thresh = thresh,
                      size = 2,
                      numSol = 1000000)
  
  subc$parameters = list.labs.association
  save(subc, file = paste(dataset$name,"_subc.RData", sep = ""))
  
  #Here the foreign classification results are saved in a .txt file 
  sink( paste(dataset$name,"_ForeignClassification_output.txt"))
  for(i in rownames(subc[["size.2"]])){
    print(confusionTable(predMap, cascade = i,other.classes='all', sort = TRUE))
  }
  sink()
    
    
  return(list(predMap,subc))
  
}



###########################################################################################
# Function for computing the confusion Matrix 
    # Input:  predMap - depicts the results of the binary classification training 
    #                   using a cross validation set up, as  well as the evaluation of all samples 
    #                   using the trained classifiers
    #         subcas - depicts the sensitivities of the pairwise binary classification, 
    #                  as well as the association between original label names and the renamed labels
    # Output: a list of all possible confusion matricies 
###########################################################################################

confusionMatrix <- function(predMap, subcas){
  labs <- predMap$meta["label",]
  
  combs <- expand.grid(0:max(labs),0:max(labs))
  combs <-combs[combs[1]<combs[2],,drop=FALSE]
  combs <- paste0("[",combs[,1],"vs",combs[,2],"]")
  
  pred <- predMap$pred
  
  # Counting the number of different label desicions, for each sample   
  predsAll <- apply(pred[combs,,drop=FALSE],2,function(x){
    sapply(0:10,function(i){sum(x==i)})
  })
  
  
  
  # checking for non-unique majorities 
  vecForNonUnqiueMaj <- list()
  sampleForNonUnqiueMaj <- c()
  idx <- 1
  for(i in 1:(dim(predsAll)[2])){
    maj <- max(predsAll[,i])
    allMaj <- which(predsAll[,i] == maj)
    
    if(length(allMaj) > 1){
      vecForNonUnqiueMaj[[idx]] <- allMaj
      sampleForNonUnqiueMaj[idx] <- i
      idx <- idx +1
    }
  }
  
  preds <- apply(predsAll,2,function(x){which.max(x)[1]-1})
  
  
  if(idx > 1){
    # If there are non-unique majorities we return all possible confusion matricies 
    
    allCombis <- expand.grid(vecForNonUnqiueMaj)
    listForAllConfMat <- list()
    
    for(k in 1:(dim(allCombis)[1])){
      preds[sampleForNonUnqiueMaj] <- as.numeric(allCombis[k,])-1
      
      confMat <- as.matrix(table(preds,labs))
      
      for(i in 1:ncol(confMat)){
        confMat[,i] <- round(confMat[,i]/sum(confMat[,i]), digits=3)
      }
      ogLabs <- sub(" = ([0-9])*","",as.vector(unlist(subcas$parameters)))
      colnames(confMat) <- ogLabs
      rownames(confMat) <- ogLabs
      
      
      
      confMat <- t(confMat)
      # confMat <-t(confMat)[c(2,3,1,4), c(4,1,3,2)] # for printing the output of confusion matrix as depicted for dataset GSE117851 
      
      samplesWithSharedMaj <- sampleForNonUnqiueMaj
      associatedLabs <- ogLabs[as.numeric(allCombis[k,])]
      
      listForAllConfMat[[k]] <- list(confMat, samplesWithSharedMaj, associatedLabs)
    }
    
    return(listForAllConfMat)
    
  }else{
    confMat <- as.matrix(table(preds,labs))
    
    for(i in 1:ncol(confMat)){
      confMat[,i] <- round(confMat[,i]/sum(confMat[,i]), digits=3)
    }
    ogLabs <- sub(" = ([0-9])*","",as.vector(unlist(subcas$parameters)))
    colnames(confMat) <- ogLabs
    rownames(confMat) <- ogLabs
    # confMat <- t(confMat)[ c(3,1,2,4), c(4,2,1,3)] # for printing the output of confusion matrix as depicted for dataset GSE73338 
    # return(confMat) # for printing the output as in the Confusion matrix for dataset GSE73338 
    return(t(confMat))
  }

}




###########################################################################################
# Computing the SVM as well as the confusion matricies for GSE73338 and GSE117851 
###########################################################################################


ovoSVM.73338 <- start.run(GSE73338)
confMat.73338 <- confusionMatrix(ovoSVM.73338[[1]], ovoSVM.73338[[2]])

# For saving the computed confusion matrix as a .txt file 
sink( "GSE73338_ConfMat_output.txt")
print(confMat.73338)
sink()





ovoSVM.117851 <- start.run(GSE117851)
confMat.117851 <- confusionMatrix(ovoSVM.117851[[1]], ovoSVM.117851[[2]])



# For saving the computed confusion matricies as a .txt file 
sink( "GSE117851_ConfMat_output.txt")
print(confMat.117851)
sink()