# Este programa recibe un fichero de prevalencia y devuelve el mismo fichero con una columna más de cluster

# cargamos el fichero de funciones auxiliares
source("/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/Henar/funciones.R")
source("/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/Henar/outliers.R")
source("/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/Henar/ordiplot.R")

# Definir numero de clusters 
num_clusters<-30
# Definir numero de iteraciones bootstraping
num_bootstrap<-10
# Definir threshold de porcentaje intra-cluster que quiero mostrar en las estadisticas finales
threshold_per_cluster <- 15

#siadh
#finput="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-phewas-prevalencia2.csv"

#hipo>18
finput="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/mayor18-total-phewas-prevalencia.csv"
mydata=read.csv(finput, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

######################### QUITAMOS LOS INDIVIDUOS QUE NO QUEREMOS ANALIZAR ###############################
# solo hombres
#mydata <- mydata[ which(mydata[4]=="Mujer"), ]

# los que son hipertensos
#mydata <- mydata[which(mydata[["401.1"]]==1),]
##########################################################################################################

# suponemos que la primera columna con codigos es la siguiente a Dx.Todos
pos_first_field<-match("Dx.Todos",names(mydata))
pos_first_field<- pos_first_field + 1

# definimos qué variables no queremos que se contemplen al calcular distancias

# en analisis de siadh podemos quitar siadh, hiponatremia, hipertension esencial, diabetes, 
# hiperlipidemia, hipercolesterolemia, tobacco disorder
codigo_quitar<-c("276.12", "401.1", "250.2", "272.11", "272.1", "318") #"253.7"

# Nos quedamos con la feature matrix quitando las columnas seleccionadas y las filas todo a cero
parcial<-subsetPhewasVars(mydata, TRUE, codigo_quitar, pos_first_field)

# Defina el path al fichero de mapeo phewas
fphewas="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/Henar/PheWAS_code_translation_v1_2-ORIGINAL.txt"

# todas las columnas (11) son de tipo character
colClasses=c(rep("character",11))
mapeo=read.csv(fphewas, header=TRUE, sep="\t", colClasses=colClasses,fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE)

# contamos el numero de pacientes
num_patients<-nrow(parcial)

mij<-as.matrix(sapply(parcial,as.numeric))

# Esta operacion quitará los outliers del data set antes de calcular las distancias
#mij<-detectaOutliers(mij)


########### Todas las distancias binarias ####################

library(ade4)

d<-list()

for (i in 1:10)
  d[[i]]<-dist.binary(mij, method=i)

################################################################

library(clustsig)
# Qué medida vamos a usar en esta ejecución
distance_measure<-d[[2]]
d2 <- function(X) ade4::dist.binary(X, method = "2")
res.single.siadh<-simprof(mij, num.expected=50, num.simulated=49,
                    method.cluster="single", method.distance=d2,
                    method.transform="identity", alpha=0.05,
                    sample.orientation="row",silent=FALSE, increment=1)

library(fpc)
library(prabclus)
library(cluster)

# evaluamos los criterios internos para el número de clusters
mejor_interno<-vector()
particiones<-list()
mejor_interno<-evaluaInternal(num_clusters,distance_measure,num_bootstrap,"ward.D2")

print(paste("Según los criterios within (imprima la variable mejor_interno), el mejor número de cluster es",which(mejor_interno==max(mejor_interno)),sep=" "))


####################### EJECUTAR HASTA AQUI Y DEFINIR DE NUEVO num_clusters ############################

num_clusters<-3

pdffile="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/hombres-mayores18-simprof100-phewas-d2-wardD.pdf"
fclusters="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/hombres-k15-phewas-d2-wardD.csv"
fstats="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/stats-cluster-hombres-k15-phewas-d2-wardD.csv"

#external_crit<-evaluaExternal(particiones[[num_clusters]]$result$partition, d, num_clusters, "ward.D")

# Qué medida vamos a usar en esta ejecución
distance_measure<-d[[2]]

# calculamos el clustering jerárquico con el método X para la matriz de distancia
hcl5<-hclust(distance_measure,method="ward.D")

#pintamos dendograma y cortamos el árbol para 3 clusters que se pintan en rojo
pdf(pdffile)	
plot(res.ward1$hclust)
rect.hclust(res.ward1$hclust, h=2.1, border="red")

# agrupamos en k clusters
groups<-cutree(res.ward1$hclust, h=2.1)

#para medir la silhoutte se usa esto.
sil<-silhouette(groups,distance_measure)
summary(sil)
dev.off()


num_clusters<-res.ward1$numgroups

salida<-as.data.frame(parcial)
#salida$cluster<-particiones[[num_clusters]]$result$partition
salida$cluster<-fillSignificant(salida,res.ward1$significantclusters)
#imprimimos la visualización del cluster
miplot<-ClusterOrdiPlot(distance_measure,salida$cluster,0)
#salida$cluster<-groups
salida$NHC<-mydata[row.names(parcial),"Nº.Historia"]

# Escribimos el fichero con los individuos por cluster con su NHC
write.table(salida,fclusters,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8")

# quitamos las dos últimas columnas que tienen datos de pacientes y del cluster
salida<-salida[,-ncol(salida)]

clusters<-list()
superan<-list()
statsclusters<-read.table(text = "",col.names = names(salida), as.is = TRUE, check.names = FALSE)

todos_superan<-vector()

for (i in 1:num_clusters) {
  
  # escogemos los datos de ese cluster como un data frame
	clusters[[i]]<-as.data.frame(salida[which(salida$cluster==i),])
  
	#transformar columnas a numerico
	for(j in c(1,1:ncol(clusters[[i]]))) {
	     clusters[[i]][,j] <- as.numeric(clusters[[i]][,j])
	}
	
	print(paste("nrows=",toString(nrow(clusters[[i]]),sep = " ")))
	
	
  #clusters[[i]]<-data.frame(salida[salida$cluster==i,],check.names = FALSE)
  #clusters[[i]]<-sapply(clusters[[i]][clusters[[i]]$cluster==i,],as.integer)
	# quito la última columna con el número de cluster
	clusters[[i]]<-clusters[[i]][,-ncol(clusters[[i]])]
	clusters[[i]]["Total",1:ncol(clusters[[i]])]<-colSums(clusters[[i]][,1:ncol(clusters[[i]])])
	clusters[[i]]["Porcentaje",1:ncol(clusters[[i]])]<-clusters[[i]]["Total",1:ncol(clusters[[i]])]/(nrow(clusters[[i]])-1)*100
	# guardo los porcentajes de cada cluster para todos los codigos phewas.
	statsclusters[i+1,1:ncol(clusters[[i]])]<-clusters[[i]]["Porcentaje",1:ncol(clusters[[i]])]
  # guardo el numero de pacientes en cada cluster
	statsclusters[i+1,ncol(clusters[[i]])+1]<-nrow(clusters[[i]])-2
	# miro qué patologías superan un threshold definido intra-cluster y guardo su posición en superan[[i]]
	superan[[i]]<-vector()
	for (j in 1:ncol(clusters[[i]])) {
	  if (statsclusters[i+1,j]>threshold_per_cluster) #guardamos su posicion
	    superan[[i]]<-c(superan[[i]],j)
	}
	# damos nombre a los elementos del vector (named num vector) qué códigos corresponden a las columnas que lo superan
	names(superan[[i]])<-colnames(parcial)[superan[[i]]]
	todos_superan<-c(todos_superan,names(superan[[i]]))
	# actualizamos las celdas remplazando las posiciones con los valores 
	for (x in 1:length(superan[[i]])) {
	  pos<-superan[[i]][x]
	  superan[[i]][x]<-statsclusters[i+1,pos]
	}
	
	# ordenamos por las patologías que tienen un porcentaje intra-cluster más alto
	superan[[i]]<-sort(superan[[i]], decreasing = T)
	
	# inicializamos la columna que albergará la cadena que describe los que superan el threshold
	statsclusters[i+1,ncol(clusters[[i]])+2]<-""
	
	# volvemos a actualizar las celdas precediéndolas con el significado del código phewas 
	for (x in 1:length(superan[[i]])) {
	  pos<-superan[[i]][x]
	  superan[[i]][x]<-paste(mapeo[which(mapeo[["phewas_code"]]==names(superan[[i]][x])),"phewas_string"][1], 
	                         as.character(round(as.numeric(superan[[i]][x]),2)), 
	                         sep = " (")
	  superan[[i]][x]<-paste(superan[[i]][x],")",sep = "")
	  # concatenamos cada significado phewas de los que superen el threshold
	  statsclusters[i+1,ncol(clusters[[i]])+2]<-paste(statsclusters[i+1,ncol(clusters[[i]])+2],superan[[i]][x], sep=" | ")
	}
}

# Rellenamos el significado de cada código phewas en la fila 1 de stats
#for (i in 1:ncol(clusters[[2]]))
#  statsclusters[1,i]<-mapeo[which(mapeo[["phewas_code"]]==colnames(statsclusters)[i]),"phewas_string"][1]             

# Escribimos el fichero con las estadísticas por cluster
write.table(statsclusters,fstats,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")




