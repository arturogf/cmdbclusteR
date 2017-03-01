#Removes everything before starting
#rm (list=ls())

library("ggplot2")

finput<-fclusters
mydata = read.csv(finput, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE)

#Number of column were de number of cluster is stored
num_clusters <- length(unique(mydata$cluster))
#Number of column were de number of GRD is stored
num_GRD <- length(unique(mydata$GRD))

#We create a table with the data for the heatmap
classes<-c("character","character","numeric")
heatdata <- read.table(text="",col.names=c("cluster","GRD","heat"), colClasses=classes)

#For each number of cluster we have the number of GRD codes and the frecuency.
for(i in 1:num_clusters){
mydataCluster <- as.data.frame(mydata[which(mydata$cluster == i),])
grd <- mydataCluster$GRD
grd.freq<- table(grd)
dtgrd<-as.data.frame(grd.freq)
dtgrd<-dtgrd[order(-dtgrd[,2]),]
old = options(digits=1)
dtgrd[,2]<-dtgrd[,2]/nrow(mydataCluster)*100

#If we have only one GRD code in one cluster we don´t generate a barplot
if(nrow(grd.freq)==1){
  print("This cluster presents only one GRD")
}else{
  # -------- define output files --------
  nombre_generado <-
    paste(
      directory, 
      "/data/processed/",
      "Cluster",
      i,
      "grafica.png",
      sep = ""
    )
  
#We obtain a png barplot for each cluster
  png(nombre_generado)
  barplot(dtgrd[1:10,2],names.arg = substr(dtgrd[1:10,1],1,20),
        main=paste("GRD más frecuentes para cluster",i),ylab="Frecuencia(%)",las=2)
  dev.off()
}


# if there is less than 10 GRD codes for this cluster
lengrd<-nrow(dtgrd)
if (lengrd>=10)
  lengrd<-10

#We fill the heatdata
a<-expand.grid(cluster=paste("cluster",i,sep = ""),GRD=unique(dtgrd[1:lengrd,1]))
print(a)
a$heat<-dtgrd[1:lengrd,2]
print(a)
heatdata<-rbind(heatdata,a)
}

# -------- define output files --------
nombre_generado_heat <-
  paste(
    directory, 
    "/data/processed/",
    "Cluster_",
    "heatmap.png",
    sep = ""
  )

png(nombre_generado_heat)

#We generate a heatmap with de prevalent GRD codes
ggplot(heatdata, aes(x = cluster, y = substr(heatdata[["GRD"]],1,20))) + 
  geom_tile(aes(fill = heat), colour="white") + 
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  xlab("Number of cluster") +
  ylab("GRD codes") +
  ggtitle("GRD codes prevalence for each cluster") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_line(colour = "white")) +
  geom_vline(xintercept = seq(0.5,num_clusters+0.5), colour = "grey") +
  geom_hline(yintercept = seq(0.5,num_GRD+0.5), colour = "grey")

dev.off()
