
`calcMeanAngles` <-
function(x, y)
{
	numSam <- length(x)
	allResults = double(numSam)
	#cat("Starting C code\n")
	print(x)
	print(y)

	print(allResults)

	result <- .C("calcMeanAnglesC",
	   as.double(x),					# X Vec 
	   as.double(y),					# Y Vec
	   as.integer(numSam),					# numSamples
	   allResults = allResults)$allResults			# results
	
	return(result)	
}



`calcAnglesToNorm` <-
function(angles)
{
	numSam <- length(angles)
	allResults = double(2*numSam)
	#cat("Starting C code\n")
	#print(x)
	#print(y)

	print(angles)

	print(allResults)

	result <- .C("calcAnglesToNormC",
	   as.double(angles),					# X Vec 
	   as.integer(numSam),					# numSamples
	   allResults = allResults)$allResults			# results
	
	result <- matrix(result, ncol = 2)

	return(result)	
}

