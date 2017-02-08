# This script is used to define all the packages used in cmbdclusteR 
# and install them previously to execution if they are not installed yet
# Note: The script should be called before executing cmbdclusteR

list.of.packages <- c("ade4","clustsig","fpc","prabclus","cluster","clusterCrit","DMwR","HighDimOut")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)