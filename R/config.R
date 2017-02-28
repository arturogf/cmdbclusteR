# This script is used to define all the packages used in cmbdclusteR 
# and install them previously to execution if they are not installed yet
# Note: The script should be called before executing cmbdclusteR
rm(list=ls())

list.of.packages <- c("ade4","clustsig","fpc","prabclus","cluster","clusterCrit","DMwR","HighDimOut", "shiny", "ggplot2","vegan")

if(as.numeric(substr(version$minor,1,1))>=3){
  list.of.packages[length(list.of.packages)+1] <- "tcltk2"
}else{
  list.of.packages[length(list.of.packages)+1] <- "tcltk"
}

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(as.numeric(substr(version$minor,1,1))>=3){
  library("tcltk2")
}else{
  library("tcltk")
}

library("shiny")

# Select directory where cmbdclusteR-master was downloaded to
print("A window will pop up to select folder where cmbdclusteR-master was downloaded to.")
directory <- tk_choose.dir(default="", caption="Select folder where cmbdclusteR-master was downloaded to")

source(file.path(directory, "/R/declaracion.R"))
