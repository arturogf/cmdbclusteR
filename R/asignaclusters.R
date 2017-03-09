# This script receives a population file with prevalence info 
# and returns the same file including a column with the cluster assigned

# include source files
source(file.path(directory, "/R/funciones.R"))
source(file.path(directory, "/R/ordiplot.R"))

# Define num of clusters
# print("Select number of clusters desired or press ENTER for default value (30):")
# num_clusters <- as.numeric(readLines(n=1, ok=FALSE))
# if (is.na(num_clusters)){
#   num_clusters <- 30
# }
# Define num of bootstraping iterations
print("Select number of bootstraping iterations desired or press ENTER for default value (500):")
num_bootstrap <- as.numeric(readLines(n=1, ok=FALSE))
if (is.na(num_bootstrap)){
  num_bootstrap <- 500
}
itfilter <- paste("simprof", num_bootstrap, "it", sep="")
# Define threshold for intra-cluster contribution that I want to show in final output
threshold_per_cluster <- 0 #REVISAR, pedir por pantalla

# Redefine mydata
mydata <- prevalentes
remove(prevalentes)

# we suppose that first column with ICD9 separated codes is the following to Dx.Todos
pos_first_field <- pos_first_ICD9
  
# -------- Define variables that we do not want to consider when calculating distances (similarity) ---------

# e.g. for SIADH, we can remove, hiponatremia, essential hipertension, diabetes, 
# hiperlipidemia, hipercolesterolemia, tobacco disorder

# Posible codes to drop in Hiponatremia  (SIADH: ICD9 = 253.6, PHEWAS=253.7)
codes_to_drop <- c("276.12", "401.1", "250.2", "272.11", "272.1", "318", "253.7")

# Posible codes to drop in mood disorder (ICD9 = 296)
#codes_to_drop <- c("401.9","272.0","250.0")
#codes_to_drop<-c("FFF","999")
#codes_to_drop<-c("")

# We subset the feature matrix, dropping the selected variables to remove and the rows with all columns==0
output <- subsetPhewasVars(mydata, TRUE, codes_to_drop, pos_first_field)
parcial <- as.data.frame(output[[1]])
mycopy <- as.data.frame(output[[2]])
remove(output)

# convert data frame to numeric matrix
mij <- as.matrix(sapply(parcial, as.numeric))

# This line would remove outliers from data before calculating distances. 
# It did not work very well with SIADH data.
#mij<-detectaOutliers(mij)

########### Calculation of Binary Distances ####################

library(ade4)

# d <- list()
# for (i in 1:10)
#   d[[i]] <- dist.binary(mij, method = i)

# Select binary distance method
cat("1 = JACCARD index (1901) S3 coefficient of GOWER & LEGENDRE\n")
cat("s1 = a/(a+b+c) --> d = sqrt(1 - s)\n")
cat("2 = SOCKAL & MICHENER index (1958) S4 coefficient of GOWER & LEGENDRE \n")
cat("s2 = (a+d)/(a+b+c+d) --> d = sqrt(1 - s)\n")
cat("3 = SOCKAL & SNEATH(1963) S5 coefficient of GOWER & LEGENDRE\n")
cat("s3 = a/(a+2(b+c)) --> d = sqrt(1 - s)\n")
cat("4 = ROGERS & TANIMOTO (1960) S6 coefficient of GOWER & LEGENDRE\n")
cat("s4 = (a+d)/(a+2(b+c)+d) --> d = sqrt(1 - s)\n")
cat("5 = CZEKANOWSKI (1913) or SORENSEN (1948) S7 coefficient of GOWER & LEGENDRE\n")
cat("s5 = 2*a/(2*a+b+c) --> d = sqrt(1 - s)\n")
cat("6 = S9 index of GOWER & LEGENDRE (1986)\n")
cat("s6 = (a-(b+c)+d)/(a+b+c+d) --> d = sqrt(1 - s)\n")
cat("7 = OCHIAI (1957) S12 coefficient of GOWER & LEGENDRE\n")
cat("s7 = a/sqrt((a+b)(a+c)) --> d = sqrt(1 - s)\n")
cat("8 = SOKAL & SNEATH (1963) S13 coefficient of GOWER & LEGENDRE\n")
cat("s8 = ad/sqrt((a+b)(a+c)(d+b)(d+c)) --> d = sqrt(1 - s)\n")
cat("9 = Phi of PEARSON = S14 coefficient of GOWER & LEGENDRE\n")
cat("s9 = ad-bc)/sqrt((a+b)(a+c)(b+d)(d+c)) --> d = sqrt(1 - s)\n")
cat("10 = S2 coefficient of GOWER & LEGENDRE\n")
cat("s10 =  a/(a+b+c+d) --> d = sqrt(1 - s) and unit self-similarity\n")
print("Select binary distance method or press ENTER for value by default (2):")
nmethod <- as.integer(readLines(n=1, ok=FALSE))
if (is.na(nmethod)){
  nmethod <- 2
}
distancefilter <- paste("d", nmethod, sep="")

distance_measure <- dist.binary(mij, method=nmethod)

################################################################

library(clustsig)
# ---------- Define which distance measure to use in this execution -------------
#distance_measure <- d[[2]]

# Identify the distance to use in simprof
d <- function(X) ade4::dist.binary(X, method = nmethod)

# --- here we have to define which clustering operations to perform, including method (ward, single, etc.) ---
# Carry out the similarity profile analysis with 50 bootstrapping iterations
res.ward.siadh <-simprof(mij, num.expected=num_bootstrap, num.simulated=(num_bootstrap-1),
                         method.cluster="ward.D", method.distance=d,
                         method.transform="identity", alpha=0.05,
                         sample.orientation="row",silent=FALSE, increment=10)

library(fpc)
library(prabclus)
library(cluster)

#create a file to save output data
if (!file.exists(file.path(directory, "data/output"))){
  dir.create(file.path(directory, "data/output"))
}

# -------- define output files --------
clusterfilter <- "ward.D" #REVISION, pedir por pantalla
nombre_generado <- paste(substring(nombre_generado, 1, nchar(nombre_generado)-4), distancefilter, clusterfilter, itfilter, sep="-")
stayboxplotfile <- file.path(directory, "data/output", paste(nombre_generado, "LOS-boxplot", ".pdf", sep=""))
pdffile <- file.path(directory, "data/output", paste(nombre_generado, ".pdf", sep=""))
ordifile <- file.path(directory, "data/output", paste(nombre_generado, "-ordiplot.pdf", sep=""))
fclusters <- file.path(directory, "data/output", paste(nombre_generado, ".csv", sep=""))
fstats <- file.path(directory, "data/output", paste(nombre_generado, "-stats.csv", sep=""))

# plot dendogram and red line for cut-off height h
pdf(pdffile)
simprof.plot(res.ward.siadh)
dev.off()

# num_clusters is taken from the output obtained from simprof execution
num_clusters <- res.ward.siadh$numgroups

salida <- as.data.frame(parcial) 

#add new column 'cluster' to salida using a function defined in funciones.R
salida$cluster <- fillSignificant(salida, res.ward.siadh$significantclusters)

#print ordiplot for cluster visualization
#pdf(ordifile)
#miplot <- ClusterOrdiPlot(distance_measure, salida$cluster, 0)
#dev.off()
#salida$cluster<-groups


#################################################################################


clusters <- list()
superan <- list()
statsclusters <- read.table(text = "",col.names = names(salida), as.is = TRUE, check.names = FALSE)

for (i in 1:num_clusters) {
  
  # select the data for cluster i as a data frame
  clusters[[i]] <- as.data.frame(salida[which(salida$cluster == i), ])
  
  # transform columns to numeric
  for (j in c(1, 1:ncol(clusters[[i]]))) {
    clusters[[i]][, j] <- as.numeric(clusters[[i]][, j])
  }
  
  # print the number of records included in each cluster
  print(paste("nrows in cluster ", toString(i), ":", toString(nrow(clusters[[i]]), sep = "")))
  
  #clusters[[i]]<-data.frame(salida[salida$cluster==i,],check.names = FALSE)
  #clusters[[i]]<-sapply(clusters[[i]][clusters[[i]]$cluster==i,],as.integer)
  
  # drop the last column with the cluster number and add the Total and Percentage rows at the end
  clusters[[i]] <- clusters[[i]][, -ncol(clusters[[i]])]
  clusters[[i]]["Total", 1:ncol(clusters[[i]])] <- colSums(clusters[[i]][, 1:ncol(clusters[[i]])])
  clusters[[i]]["Porcentaje", 1:ncol(clusters[[i]])] <- clusters[[i]]["Total", 1:ncol(clusters[[i]])] / (nrow(clusters[[i]]) - 1) *100
  
  # store the percentages for each cluster for all the phewas codes
  statsclusters[i + 1, 1:ncol(clusters[[i]])] <-0
  statsclusters[i + 1, 1:ncol(clusters[[i]])] <-clusters[[i]][nrow(clusters[[i]]), 1:ncol(clusters[[i]])]
  
  # store the number of episodes in each cluster 
  statsclusters[i + 1, ncol(clusters[[i]]) + 1] <- nrow(clusters[[i]]) - 2
  
  # We create a new column of 0 to store de number of unique patients in each cluster
  statsclusters[i + 1, ncol(clusters[[i]]) + 2] <- 0
  
  # look which pathologies' presence is over intra-cluster threshold and store its position in superan[[i]]
  superan[[i]] <- vector()
  for (j in 1:ncol(clusters[[i]])) {
    if (statsclusters[i + 1, j] > threshold_per_cluster)
      superan[[i]] <- c(superan[[i]], j)
  }
  # give names to superan[[i]] vector elements using the codes from original data 
  names(superan[[i]]) <- colnames(parcial)[superan[[i]]]
  
  # update cells of superan[[i]] by replacing stored positions with values 
  for (x in 1:length(superan[[i]])) {
    pos <- superan[[i]][x]
    superan[[i]][x] <- statsclusters[i + 1, as.numeric(pos)]
  }
  
  # order superan[[i]] in decreasing order of intra-cluster contribution
  superan[[i]] <- sort(superan[[i]], decreasing = T)
  
  # initialize column that will store description for the codes that are higher than defined threshold
  statsclusters[i + 1, ncol(clusters[[i]]) + 3] <- ""
  
  if (input_type=="PHEWAS"){
    # we update superan[i] cells again, preceeding them with each phewas code description
    for (x in 1:length(superan[[i]])) {
      pos <- superan[[i]][x]
      superan[[i]][x] <-
        paste(mapeo[which(mapeo[["phewas_code"]] == names(superan[[i]][x])), "phewas_string"][1],
              as.character(round(as.numeric(superan[[i]][x]), 2)),
              sep = " (")
      superan[[i]][x] <- paste(superan[[i]][x], ")", sep = "")
      
      # concatenate each phewas description for all codes that are higher than threshold 
      statsclusters[i+1,ncol(clusters[[i]])+3]<-paste(statsclusters[i+1,ncol(clusters[[i]])+3],superan[[i]][x], sep=" | ")
    }
    for (k in 1:(length(statsclusters)-3)){
      posi <- which(mapeo[,3] == (names(statsclusters[,1:(length(statsclusters)-3)])[k]))
      if (length(posi))
        statsclusters[1,k]<- mapeo[posi[1],4]
      else
        statsclusters[1,k]<- names(statsclusters[,1:(length(statsclusters)-3)])[k]
    }
  }
  else {
    # we update superan[i] cells again, preceeding them with each phewas code description
    for (x in 1:length(superan[[i]])) {
      pos <- superan[[i]][x]
      superan[[i]][x] <-
        paste(mapeo[which(mapeo[["icd9"]] == names(superan[[i]][x])), "icd9_string"][1],
              as.character(round(as.numeric(superan[[i]][x]), 2)),
              sep = " (")
      superan[[i]][x] <- paste(superan[[i]][x], ")", sep = "")
      
      # concatenate each phewas description for all codes that are higher than threshold 
      statsclusters[i+1,ncol(clusters[[i]])+3]<-paste(statsclusters[i+1,ncol(clusters[[i]])+3],superan[[i]][x], sep=" | ")
    }
    # rellenar con etiquetas de icd9_string la segunda fila
    for (k in 1:(length(statsclusters)-3)){
      posi <- which(mapeo[,1] == (names(statsclusters[,1:(length(statsclusters)-3)])[k]))
      if (length(posi))
        statsclusters[1,k]<- mapeo[posi,2]
      else
        statsclusters[1,k]<-names(statsclusters[,1:(length(statsclusters)-3)])[k]
    }
  }
}



# Rellenamos el significado de cada código phewas en la fila 1 de stats
#for (i in 1:ncol(clusters[[2]]))
#  statsclusters[1,i]<-mapeo[which(mapeo[["phewas_code"]]==colnames(statsclusters)[i]),"phewas_string"][1] 

# ---- if needed, add the patient record number if needed for exploration -----
#salida$NHC <- mydata[row.names(parcial), pos_numeroHC]

colnames(statsclusters)[ncol(statsclusters)] <- "Intracluster_contribution"
colnames(statsclusters)[ncol(statsclusters)-1] <- "unique_patients"
colnames(statsclusters)[ncol(statsclusters)-2] <- "nepisodes"

# ---- if needed, add the GRD number for exploration in salida -----
#salida$GRD <- mydata[row.names(parcial), pos_GRD]

salida <- cbind(mycopy, salida)

# we create labels for each cluster where we include in brackets the number of episodes
labls <- c()
for (i in seq(1:length(unique(salida$cluster))))
  labls <-
  c(labls,
    paste(
      "cluster",
      as.character(i),
      "[n=",
      as.character(length(which(salida$cluster == i))),
      "]",
      sep = ""
    ))

library("ggplot2")
#diagnostico1<-factor(substr(salida$D1, 1, 3))

pdf(stayboxplotfile, paper = "a4r")

#plot the boxplot
ggplot(salida,
       aes(
         group = salida$cluster,
         x = salida$cluster,
         y = salida[, pos_stay]
       ))  +
  coord_flip() +
  geom_boxplot(
   #aes(fill = diagnostico1),
    #fill="white",
    colour = "black",
    outlier.colour = "red",
    outlier.shape = 2
  ) +
  stat_boxplot(geom = 'errorbar') +
  geom_jitter(width = 0.2, height=0.2) +
  #geom_point(position=position_dodge(width=0.75),aes(group=salida$cluster))+
  labs(title = "Boxplot: LOS per cluster") +
  ylab("Length of Stay(LOS)") +
  scale_x_continuous(name = "Num. cluster",
                     breaks = seq(1:length(unique(salida$cluster))),
                     labels = labls)

dev.off()

# We write the file with all the statitics per cluster
write.table(salida,fclusters,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")


# We write the file with all the statitics per cluster
write.table(salida,fclusters,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")

for (i in 1:num_clusters) {
  # store the number of unique patients in each cluster
  statsclusters[i + 1, ncol(clusters[[i]]) + 2] <- length(unique(as.data.frame(salida[which(salida$cluster == i), ])[,pos_numeroHC]))
  #print(paste("i: ",statsclusters[i + 1, ncol(clusters[[i]]) + 2]))
}

# we create two different dataframes to be able to convert the 2nd to numeric        
statsclusters1<-statsclusters[1,]
statsclusters2<-statsclusters[2:nrow(statsclusters),]
statsclusters2[,1:(ncol(statsclusters2)-3)]<-sapply(statsclusters2[,1:(ncol(statsclusters2)-3)], as.numeric)

# We write the file with all the statitics per cluster
write.table(statsclusters1,fstats,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")
write.table(statsclusters2,fstats,append=TRUE,sep=";",row.names = FALSE, col.names=FALSE, fileEncoding = "UTF-8",dec=",")
