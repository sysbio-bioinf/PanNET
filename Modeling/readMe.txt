Instructions for running the modeling part:

The script will require you to set up a working directory (setwd("...")) and requires dependencies from a set of working packages listed in the first lines of the code.
Each subpart of the script is divided by the figure it corresponds to from either main or supplement.
When loading the model, it is possible to choice between the delayed and non-delayed network. In the delayed model, some nodes providing only linear connections are represented as time delays. This helps to reduce the computational requirements of some simulations. Both models will return the same patterns of attractors.
There are some analyses that require higher computational power (indicated in the code). For this reason, we also provided the resulting RData-files to allow a quick reproduction of the results. This aspect affects the tumor driver analysis and the basin evaluation.