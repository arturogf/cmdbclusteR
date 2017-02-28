#Removes everything before starting
rm (list=ls())

# Select directory where cmbdclusteR-master was downloaded to
print("A window will pop up to select folder where cmbdclusteR-master was downloaded to.")
directory <- tk_choose.dir(default="", caption="Select folder where cmbdclusteR-master was downloaded to")

#input
finput="C:/Users/user/Documents/2016-2017/Practicas_Clinico/Codigos/Output/simprof1000-phewas-d2-wardD_mood_women.csv"
mydata = read.csv(finput, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE)

num_clusters <- length(unique(mydata$cluster))

classes<-c("character","character","numeric")
heatdata <- read.table(text="",col.names=c("cluster","GRD","heat"), colClasses=classes)

for(i in 1:num_clusters){
mydataCluster <- as.data.frame(mydata[which(mydata$cluster == i),])
grd <- mydataCluster$GRD
grd.freq<- table(grd)
#grd.freq <- grd.freq[order(grd.freq[,2])]
dtgrd<-as.data.frame(grd.freq)
dtgrd<-dtgrd[order(-dtgrd[,2]),]
old = options(digits=1)
dtgrd[,2]<-dtgrd[,2]/nrow(mydataCluster)*100

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

  png(nombre_generado)
  barplot(dtgrd[1:10,2],names.arg = substr(dtgrd[1:10,1],1,20),
        main=paste("GRD más frecuentes para cluster",i),ylab="Frecuencia(%)",las=2)
  dev.off()
}
a<-expand.grid(cluster=paste("cluster",i,sep = ""),GRD=unique(dtgrd[,1]))
a$heat<-dtgrd[,2]
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
ggplot(heatdata, aes(x = cluster, y = substr(heatdata[["GRD"]],1,20))) + 
  geom_tile(aes(fill = heat), colour="white") + 
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  xlab("Number of cluster") +
  ylab("GRD codes") +
  ggtitle("GRD codes prevalence for each cluster") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_line(colour = "grey")) +
  geom_vline(xintercept = seq(0.5,num_clusters+0.5), colour = "red")

  #Hacer las horizontales también y quitar las negras

dev.off()
