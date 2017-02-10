# This script receives a population file with prevalence info 
# and returns the same file including a column with the cluster assigned

# include source files
source("./funciones.R")
#source("./outliers.R")
source("./ordiplot.R")

# Define num of clusters
num_clusters <- 30
# Define num of bootstraping iterations
num_bootstrap <- 10
# Define threshold for intra-cluster contribution that I want to show in final output
threshold_per_cluster <- 15

#siadh with prevalence information
finput="/Users/arturogf/cmdbclusteR/data/processed/siadh-phewas-prevalencia2.csv"
mydata = read.csv(finput, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

######################### Remove individuals that we do not want to analyze ##############################
# only women
#mydata <- mydata[ which(mydata[4]=="Mujer"), ]

# only hypertensive
#mydata <- mydata[which(mydata[["401.1"]]==1),]
##########################################################################################################

# we suppose that first column with ICD9 separated codes is the following to Dx.Todos
pos_first_field <- match("Dx.Todos", names(mydata))
pos_first_field<- pos_first_field + 1

# -------- Define variables that we do not want to consider when calculating distances (similarity) ---------

# e.g. for SIADH, we can remove, hiponatremia, essential hipertension, diabetes, 
# hiperlipidemia, hipercolesterolemia, tobacco disorder
codes_to_drop <- c("276.12", "401.1", "250.2", "272.11", "272.1", "318") #"253.7"

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
res.ward.siadh <-simprof(mij, num.expected=500, num.simulated=499,
                    method.cluster="ward", method.distance=d2,
                    method.transform="identity", alpha=0.05,
                    sample.orientation="row",silent=FALSE, increment=1)

library(fpc)
library(prabclus)
library(cluster)

# With this lines we can evaluate internal ('within') criteria for clusters
best_within <- vector()
particiones <- list()
best_within <-
  evaluaInternal(num_clusters, distance_measure, num_bootstrap, "ward.D2")

print(
  paste(
    "Según los criterios within (imprima la variable best_within), el mejor número de cluster es",
    which(best_within == max(best_within)),
    sep = " "
  )
)

####################### The following piece of code was used to compare output with hclust method: needs refactoring! ########################
####################### Here we stop execution and define again num_clusters, after deciding which num. is best ############################

num_clusters<-3

pdffile="/Users/arturogf/cmdbclusteR/data/output/hombres-mayores18-simprof100-phewas-d2-wardD.pdf"
fclusters="/Users/arturogf/cmdbclusteR/data/output/hombres-k15-phewas-d2-wardD.csv"
fstats="/Users/arturogf/cmdbclusteR/data/output/stats-cluster-hombres-k15-phewas-d2-wardD.csv"

#external_crit<-evaluaExternal(particiones[[num_clusters]]$result$partition, d, num_clusters, "ward.D")

# Qué medida vamos a usar en esta ejecución
distance_measure<-d[[2]]

# calculate hierarchical clustering using ward method and hclust function
hcl5 <- hclust(distance_measure, method = "ward.D")

# plot dendogram and red line for cut-off height h
pdf(pdffile)
plot(res.ward.siadh$hclust)
rect.hclust(res.ward.siadh$hclust, h = 2.1, border = "red")

# group cutting tree with height h
groups <- cutree(res.ward.siadh$hclust, h = 2.1)

# measuring silhouette of groups
sil <- silhouette(groups, distance_measure)
summary(sil)
dev.off()

############################################ end refactoring ################################################

# num_clusters is taken from the output obtained from simprof execution
num_clusters <- res.ward.siadh$numgroups

salida <- as.data.frame(parcial)
#salida$cluster<-particiones[[num_clusters]]$result$partition

#add new column 'cluster' to salida using a function defined in funciones.R
salida$cluster <- fillSignificant(salida, res.ward.siadh$significantclusters)

#print ordiplot for cluster visualization
miplot <- ClusterOrdiPlot(distance_measure, salida$cluster, 0)

#salida$cluster<-groups

# ---- if needed, add the patient record number if needed for exploration -----
salida$NHC <- mydata[row.names(parcial), "Nº.Historia"]

# We write the file with individuals and cluster assignment plus patient record number
write.table(salida,fclusters,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8")

#################################################################################

# we remove the last two columns including the cluster assignment and patient record number
salida <- salida[, -ncol(salida)]

clusters <- list()
superan <- list()
statsclusters <- read.table(text = "",col.names = names(salida), as.is = TRUE, check.names = FALSE)

todos_superan <- vector()

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
	
	# store the number of patients in each cluster
	statsclusters[i + 1, ncol(clusters[[i]]) + 1] <- nrow(clusters[[i]]) - 2
	
	# look which pathologies' presence is over intra-cluster threshold and store its position in superan[[i]]
	superan[[i]] <- vector()
	for (j in 1:ncol(clusters[[i]])) {
	  if (statsclusters[i + 1, j] > threshold_per_cluster)
	    superan[[i]] <- c(superan[[i]], j)
	}
	# give names to superan[[i]] vector elements using the codes from original data and concatenate in vector todos_superan those code names
	names(superan[[i]]) <- colnames(parcial)[superan[[i]]]
	todos_superan <- c(todos_superan, names(superan[[i]]))
	
	# update cells of superan[[i]] by replacing stored positions with values 
	for (x in 1:length(superan[[i]])) {
	  pos <- superan[[i]][x]
	  superan[[i]][x] <- statsclusters[i + 1, pos]
	}
	
	# order superan[[i]] co in decreasing order of intra-cluster contribution
	superan[[i]] <- sort(superan[[i]], decreasing = T)
	
	# initialize column that will store description for the codes that are higher than defined threshold
	statsclusters[i + 1, ncol(clusters[[i]]) + 2] <- ""
	
	# we update superan[i] cells again, preceeding them with each phewas code description
	for (x in 1:length(superan[[i]])) {
	  pos <- superan[[i]][x]
	  superan[[i]][x] <-
	    paste(mapeo[which(mapeo[["phewas_code"]] == names(superan[[i]][x])), "phewas_string"][1],
	          as.character(round(as.numeric(superan[[i]][x]), 2)),
	          sep = " (")
	  superan[[i]][x] <- paste(superan[[i]][x], ")", sep = "")
	  
	  # concatenate each phewas description for all codes that are higher than threshold 
	  statsclusters[i+1,ncol(clusters[[i]])+2]<-paste(statsclusters[i+1,ncol(clusters[[i]])+2],superan[[i]][x], sep=" | ")
	}
}

# Rellenamos el significado de cada código phewas en la fila 1 de stats
#for (i in 1:ncol(clusters[[2]]))
#  statsclusters[1,i]<-mapeo[which(mapeo[["phewas_code"]]==colnames(statsclusters)[i]),"phewas_string"][1]             

# We write the file with all the statitics per cluster
write.table(statsclusters,fstats,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")




