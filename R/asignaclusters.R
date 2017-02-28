# This script receives a population file with prevalence info 
# and returns the same file including a column with the cluster assigned

# include source files
source(file.path(directory, "/R/funciones.R"))
source(file.path(directory, "/R/ordiplot.R"))

# Define num of clusters
num_clusters <- 30
# Define num of bootstraping iterations
num_bootstrap <- 10
# Define threshold for intra-cluster contribution that I want to show in final output
threshold_per_cluster <- 0

#siadh with prevalence information
finput="/Users/arturogf/cmdbclusteR/data/processed/siadh-phewas-prevalencia2.csv"
mydata = read.csv(finput, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

# we suppose that first column with ICD9 separated codes is the following to Dx.Todos
pos_first_field <- match("Dx.Todos", names(mydata))+1

# -------- Define variables that we do not want to consider when calculating distances (similarity) ---------

# e.g. for SIADH, we can remove, hiponatremia, essential hipertension, diabetes, 
# hiperlipidemia, hipercolesterolemia, tobacco disorder
# Posible codes to drop in Hiponatremia  (SIADH: ICD9 = 253.6, PHEWAS=253.7)
codes_to_drop <- c("276.12", "401.1", "250.2", "272.11", "272.1", "318", "253.7")
# Posible codes to drop in mood disorder (ICD9 = 296)
#codes_to_drop <- c("401.9","272.0","250.0")
#codes_to_drop<-""

# We subset the feature matrix, dropping the selected variables to remove and the rows with all columns==0
parcial <- subsetPhewasVars(mydata, TRUE, codes_to_drop, pos_first_field)

# -------- define path to PheWAS mapping file ----------
fphewas = "/Users/arturogf/cmdbclusteR/data/mappings/PheWAS_code_translation_v1_2-ORIGINAL.txt"

# all columns (11) are defined of type character
colClasses = c(rep("character", 11))
mapeo=read.csv(fphewas, header=TRUE, sep="\t", colClasses=colClasses,fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE)

# we count the number of patients
num_patients <- nrow(parcial)

# convert data frame to numeric matrix
mij <- as.matrix(sapply(parcial, as.numeric))

# This line would remove outliers from data before calculating distances. 
# It did not work very well with SIADH data.
#mij<-detectaOutliers(mij)


########### Calculation of Binary Distances ####################

library(ade4)

d <- list()

for (i in 1:10)
  d[[i]] <- dist.binary(mij, method = i)

################################################################

library(clustsig)
# ---------- Define which distance measure to use in this execution -------------
distance_measure <- d[[2]]

# Identify the distance to use in simprof
d2 <- function(X) ade4::dist.binary(X, method = "2")

# --- here we have to define which clustering operations to perform, including method (ward, single, etc.) ---
# Carry out the similarity profile analysis with 50 bootstrapping iterations
res.ward.siadh <-simprof(mij, num.expected=10, num.simulated=9,
                         method.cluster="ward.D", method.distance=d2,
                         method.transform="identity", alpha=0.05,
                         sample.orientation="row",silent=FALSE, increment=10)

library(fpc)
library(prabclus)
library(cluster)

# -------- define output files --------
pdffile="/Users/arturogf/cmdbclusteR/data/output/mayores18-simprof500-phewas-d2-wardD.pdf"
ordifile="/Users/arturogf/cmdbclusteR/data/output/ordiplot-mayores18-simprof500-phewas-d2-wardD.pdf"
fclusters="/Users/arturogf/cmdbclusteR/data/output/simprof500-phewas-d2-wardD.csv"
fstats="/Users/arturogf/cmdbclusteR/data/output/stats-simprof500-cluster-phewas-d2-wardD.csv"

# plot dendogram and red line for cut-off height h
pdf(pdffile)
simprof.plot(res.ward.siadh)
dev.off()

# num_clusters is taken from the output obtained from simprof execution
num_clusters <- res.ward.siadh$numgroups

salida <- as.data.frame(parcial)
#salida$cluster<-particiones[[num_clusters]]$result$partition

#add new column 'cluster' to salida using a function defined in funciones.R
salida$cluster <- fillSignificant(salida, res.ward.siadh$significantclusters)

#print ordiplot for cluster visualization
pdf(ordifile)
miplot <- ClusterOrdiPlot(distance_measure, salida$cluster, 0)
dev.off()
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
  
  # print the number of records
  print(paste("nrows=", toString(nrow(clusters[[i]]), sep = " ")))
  
  #clusters[[i]]<-data.frame(salida[salida$cluster==i,],check.names = FALSE)
  #clusters[[i]]<-sapply(clusters[[i]][clusters[[i]]$cluster==i,],as.integer)
  
  # drop the last column with the cluster number and add the Total and Percentage rows at the end
  clusters[[i]] <- clusters[[i]][, -ncol(clusters[[i]])]
  clusters[[i]]["Total", 1:ncol(clusters[[i]])] <- colSums(clusters[[i]][, 1:ncol(clusters[[i]])])
  clusters[[i]]["Porcentaje", 1:ncol(clusters[[i]])] <- clusters[[i]]["Total", 1:ncol(clusters[[i]])] / (nrow(clusters[[i]]) - 1) *100
  
  # store the percentages for each cluster for all the phewas codes
  statsclusters[i + 1, 1:ncol(clusters[[i]])] <-clusters[[i]]["Porcentaje", 1:ncol(clusters[[i]])]
  
  # store the number of episodes in each cluster
  statsclusters[i + 1, ncol(clusters[[i]]) + 1] <- nrow(clusters[[i]]) - 2
  
  # We create a new colon of 0 to store de numbre of unique patients in each cluster
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
    superan[[i]][x] <- statsclusters[i + 1, pos]
  }
  
  # order superan[[i]] co in decreasing order of intra-cluster contribution
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
  }
}

# Rellenamos el significado de cada código phewas en la fila 1 de stats
#for (i in 1:ncol(clusters[[2]]))
#  statsclusters[1,i]<-mapeo[which(mapeo[["phewas_code"]]==colnames(statsclusters)[i]),"phewas_string"][1] 

# ---- if needed, add the patient record number if needed for exploration -----
salida$NHC <- mydata[row.names(parcial), pos_numeroHC]

colnames(statsclusters)[ncol(statsclusters)] <- "Intracluster_contribution"
colnames(statsclusters)[ncol(statsclusters)-1] <- "unique_patients"
colnames(statsclusters)[ncol(statsclusters)-2] <- "nepisodes"

# ---- if needed, add the GRD number for exploration in salida -----
salida$GRD <- mydata[row.names(parcial), pos_GRD]

for (i in 1:num_clusters) {
  # store the number of unique patients in each cluster
  statsclusters[i + 1, ncol(clusters[[i]]) + 2] <- length(unique(as.data.frame(salida[which(salida$cluster == i), ])$NHC))
  #print(paste("i: ",statsclusters[i + 1, ncol(clusters[[i]]) + 2]))
}
# we remove the patient record number from last column
#salida <- salida[, -ncol(salida)]            

# We write the file with all the statitics per cluster
write.table(statsclusters,fstats,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")