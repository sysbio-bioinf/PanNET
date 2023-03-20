
`THRESHOLDerrorMinimizer` <-
function(data, labs, weights, direction = "b")
{
	#cat("Original labs:", labs, "\n")
	#cat("Original weights:", weights, "\n")
    
    geq <- 1
    leq <- 1
    if(direction == "l")
    {
        geq <- 0
    }
    
    if(direction == "g")
    {
        leq <- 0
    }

	numFea <- nrow(data)  #Number of Samples 

	posIndex <- which(labs==1)
	negIndex <- which(labs==0)

	num_pos	<- length(posIndex)
	num_neg <- length(negIndex)

	#cat("Sorting data\n")
	data <- data[, c(posIndex, negIndex), drop = FALSE]

	#cat("Sorting weights\n")
	sorted.weights <- weights[c(posIndex, negIndex)]

	#print(sorted.weights)

	sorted.feature <- double(ncol(data))
	
	#cat("generating normVecs\n")
	thresholds <- double(ncol(data)+1)
	
	#print(normVecs[1:10,])

	allResults = double(4 * numFea)
	#cat("Starting C code\n")

	result <- .C("THRESHOLDerrorMinimizer",
	   as.double(t(data)),					#Dataset (sorted) 
	   as.double(sorted.weights),			#Weigths (sorted)
	   as.double(sorted.feature),			#Feature (sorted)
	   as.double(thresholds),				#List of thresholds
	   as.integer(numFea),					#Number of Features
	   as.integer(num_pos),					#Index of last positive example
	   as.integer(ncol(data)),				#Index of last negative example
	   as.integer(length(thresholds)),		#Index of last thresholds
       as.integer(geq),                     #Allow >= threshold (Yes = 1, No = Else)
       as.integer(leq),                     #Allow <= threshold (Yes = 1, No = Else)
	   allResults = allResults)$allResults  #The list of all results
	
	#cat("Stopping C code\n")
	
	result <- t(matrix(result,nrow = 4))
	colnames(result) <- c("Fea", "Dir", "Thr","Score")

	result[,1] <- result[,1]+1
	

	return(list(
			base.learners = result,
			weights	      = weights,
			type	      = "minimize"
	))	
}


`baselearner.threshold` <- function(data, labs,type = "error", weights = NULL, direction = "b")
{
	if(type == "error")
	{
		numSam <- length(labs)
		return(THRESHOLDerrorMinimizer(data, labs, rep(1/numSam,numSam), direction))
	} 
	if(type == "classwise")
	{
		indexPos <- labs==1
		indexNeg <- labs==0

		weights <- rep(0, length(labs))
		weights[indexPos] <- 1/(sum(indexPos)*2)
		weights[indexNeg] <- 1/(sum(indexNeg)*2)

		return(THRESHOLDerrorMinimizer(data, labs, weights, direction))
	}
	if(type == "weighted")
	{
		if(is.null(weights))
		    stop("Error in function baselearner.threshold: Argument weights is empty.")

		weights =weights/sum(weights)
		return(THRESHOLDerrorMinimizer(data, labs, weights, direction))
	}

        stop("Error in function baselearner.threshold: Unknown type selected.")
}



`prediction.threshold` <- function(models,newdata)
{
	numFea <- nrow(newdata)
	numSam <- ncol(newdata)

	features   <- models[,"Fea"]-1
	directions <- models[,"Dir"]
	thresholds <- models[,"Thr"]

	numMod <- nrow(models)

	if(max(features)+1 > numFea)
		stop("Trained classification models do not fit to newdata: Unknown feature dimension selected.")

	allResults <- rep(0, numMod*numSam )

	result <- .C("THRESHOLDprediction",
	   as.double(t(newdata)),				#Dataset (sorted) 
	   as.integer(features),				#Model parameter: Feature
	   as.integer(directions),				#Model parameter: Direction
	   as.double(thresholds),				#Model parameter: Threshold
	   as.integer(numFea),					#Number of Features (newdata)
	   as.integer(numSam),					#Number of Samples (newdata)
	   as.integer(numMod),					#Number of Modles (models)
	   allResults = as.double(allResults))$allResults	#The list of all results    

        result <- t(matrix(result, nrow = numSam, ncol = numMod ))

	return(result)
}

