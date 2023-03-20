For performing the foreign classification, as well as the multiclass Classification, in R, the following procedure was followed:

1) The local package "THRESHOLDerrorMinimizer", provided in the folder under the same name, was installed in the following way
	locallib <- tempfile()
	dir.create(locallib)
	.libPaths(locallib)
	install.packages("Path/To/THRESHOLDerrorMinimizer", lib = locallib, type = "source", repos = NULL)

2) The file linearOneVsOneSVM.R was run. 
	Here, two main functions are implemented for performing the binary classification via a linear SVM, and for constructing the multiclass desicions. 
	These procedures utilize the TunePareto, as well as the ORION package (versions of the used packages are documented in the .R file).
	Mutliple interim results are saved in either .RData files or .txt files. 
	For computing the binary SVMs the function start.run() returns a list comprising the results of the pairwise training using a crossvalidation set up and evaluation of all samples using the trained classifiers. 
	Here, before each training the selection of the 200 features with the smallest error rate considering a found threshold on the feature, is performed.  
	For constructing the confusion matricies, the function confusionMatrix(), takes the list elements returned by start.run(), as the first and second input argument, and returns all possible confusion matricies.  
	The results of the foreign classification is saved as Dataset_ForeignClassification_output.txt file, but can also be retrieved from the output of the start.run() function. 
	The results of the multiclass reclassification is returned by the confusionMatrix() function.
