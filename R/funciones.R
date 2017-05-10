#' This function subset and filter the columns indicated as well as rows with zero values
#' @param mydata this is a dataframe containing a cmbd file plus feature matrix of phewas
#' @param nozero this indicates if I want to remove rows with no patology at all
#' @param quitar_codigos a vector of variables not to be considered in the distance function
#' @param post_first_field the column that has the first Phewas code
#' @return parcial the filtered dataframe
subsetPhewasVars<-function(mydata,nozero=TRUE,dupli=FALSE,quitar_codigos, pos_first_field) {
  # seleccionamos el subconjunto de las variables icd9
  parcial<-mydata[1:(nrow(mydata)-3),pos_first_field:ncol(mydata)]
  
  # quitamos las variables que sean
  if (any(quitar_codigos!="")){
    pos_to_drop <- match(quitar_codigos,names(parcial))
    pos_to_drop <- setdiff(pos_to_drop, NA)
    if(length(pos_to_drop)!=0){
      parcial<-parcial[,-pos_to_drop]
    }
  }
  if (nozero) {
    # buscamos todas las filas que tengan todo a cero y las eliminamos. Esto permite que no haya NAs luego
    zeros<-which(apply(parcial[,], MARGIN = 1, function(x) any(x != 0))==FALSE)
    if (length(zeros)) {
      print(paste(as.character(length(zeros)),"samples were dropped before clustering because all values are abscent (zero)",sep=" "))
      parcial<-parcial[-zeros,]
      mycopy<-mycopy[-zeros,]
    }
  }
  
  if (dupli) {
    duplis<-which(duplicated(parcial)==TRUE)
    if (length(duplis)) {
      print(paste(as.character(length(duplis)),"samples were dropped before clustering because they were duplicated (all column values equal to a previous sample)",sep=" "))
      parcial<-parcial[-duplis,]
      mycopy<-mycopy[-duplis,]
    }
  }
  
  output <- list(parcial, mycopy)
  return(output)
}

#' This function returns a vector indicating on each position (for num_clusters) how many
#' criteria returned it was the best among all the numbers of clusters tried.
#' @param num_clusters the maximum number of clusters to evaluate
#' @param d the distance matrix
#' @param num_bootstrap the number of bootstrapping iterations to do (with the subset method)
#' @return mejor_interno the vector with the evaluations
evaluaInternal<-function(num_clusters,d,num_bootstrap,submethod) {

  library(clusterCrit)
  
  # Hacemos bootstrapping con subset del 50% para cada k y exploramos los criterios de assessment de clusters de la librería clusterCrit
  partitions<-list()
  criteria<-list()
  for (i in 2:num_clusters) {
    partitions[[i]] <- clusterboot(d, 
                                   B=num_bootstrap,
                                   bootmethod = "boot", 
                                   #subtuning=floor(nrow(mij)*0.2), 
                                   bscompare = TRUE,
                                   multipleboot = TRUE, 
                                   distances=TRUE,
                                   clustermethod = hclustCBI, 
                                   method=submethod, 
                                   k=i) 
    criteria[[i]]<-intCriteria(mij,partitions[[i]]$partition,"all")
  }
  
  # lo de que sea indice 2 es porque es seguro que al menos se harán 2 clusters
  internal_crit<-read.table(text = "",col.names = names(criteria[[2]]), as.is = TRUE, check.names = FALSE)
  
  for (i in 2:num_clusters) {
    internal_crit[i,]<-criteria[[i]]
  }
  
  # Miramos que num clusters es mejor
  best<-vector()
  
  for (j in 1:ncol(internal_crit)) {
    best<-c(best,bestCriterion(internal_crit[2:num_clusters,j], names(internal_crit)[j]))
  }
  
  # el mejor definitivo es un vector de conteo de cada k
  mejor_interno<-vector()
  # el primero lo ponemos a 0, ya que nunca tendremos k=1
  mejor_interno<-0
  
  # con esto contamos el número de apariciones en best. 
  #Hay que poner k-1 ya que en best se guardaban 1's para k=2, etc.
  for (k in 2:num_clusters)
    mejor_interno<-c(mejor_interno,length(which(best==k-1)))

particiones<<-partitions

return(mejor_interno)

}


evaluaExternal<-function(refpart,d, num_clusters, methodHC) {

  part<-list()
  salida<-list()
  
  for (i in 1:10) {
    # calculamos el clustering jerárquico con el método X para la matriz de distancia
    part[[i]]<-hclust(d[[i]],method=methodHC)
    
    # agrupamos en k clusters
    groups<-cutree(part[[i]], k=num_clusters)
    
    salida[[i]]<-as.data.frame(parcial)
    salida[[i]]$cluster<-groups
  }
  
  external_crit<-list()
  
  for (i in 1:10) {
        external_crit[[i]]<-extCriteria(salida[[i]]$cluster,refpart,"all")
  }

  return(external_crit)
}

#' This function takes the salida dataframe and appends the output of significant clusters
#' obtained by a simprof analysis
#' @param salida the output matrix where clustered instances are located
#' @param the lists of arrays that have the cluster found significants and their members
#' @return the column with cluster assignment for each member, ready to be added to the salida dataframe
fillSignificant<-function(salida,significantclusters) {
  
  salida$significant<-0
  
  for (i in 1:nrow(salida)) {
    
    for (j in 1:length(significantclusters))
    {
      pos<-match(i,significantclusters[[j]])
      if (!is.na(pos)) {
        salida[i,"significant"]<-j
      }
    }
  }
  
  return (salida$significant)
}

# Comparamos los clusters obtenidos con network_analysis con los obtenidos por hierarchical clustering
## -----pcor-----
compare_pcor_agnes<-function(superan,stats) {

    k<-NULL
  for (i in 1:length(superan))
    k<-c(k,paste("C_h",i,sep=""))
  
  ddd<-data.frame(k)
  
  # Lista de valores para los 11 clusters con alpha=0.05
  # cr<-list()
  # cr[[1]]<-c("800.1", "285.1")
  # cr[[2]]<-c("433.31", "571.5", "1010", "174.11")
  # cr[[3]]<-c("743.11", "743.21")
  # cr[[4]]<-c("261.2", "260", "276.5", "426.91")
  # cr[[5]]<-c("V45.77", "V88.01")
  # cr[[6]]<-c("41.9", "41.4", "591", "41")
  # cr[[7]]<-c("994.2", "509.1", "509.8", "496.21", "519.8", "600", "599.2")
  # cr[[8]]<-c("198.4", "198.6", "165.1", "198.5", "198.1", "198.2", "153.2")
  # cr[[9]]<-c("286.2", "427.21", "427.3", "979", "428.1", "585.1", "585.3", "401.22")
  # cr[[10]]<-c("411.2", "411.4")
  # cr[[11]]<-c("318", "317")
  
 
  ## Graph_pcor
   library(igraph)
  
  ig<-as.igraph(Graph_pcor, attributes = TRUE)
  cl<-clusters(ig)
  pos<-which(cl$csize>1)
  
  clus<-list()
  cr<-list()
  
  for (i in 1:length(pos)) {
    clus[[i]]<-which(cl$membership==pos[i])
    cr[[i]]<-V(ig)$label[clus[[i]]]
  }
  
 
    for (k in 1:length(cr)) {
    for (i in 1:length(superan)) {
      print(paste("cluster",i,sep=" "))
      contrib<-""
      for (j in 1:length(cr[[k]])) {
        contrib<-paste(contrib,cr[[k]][j],"(",round(superan[[i]][cr[[k]][j]],1),"), ",sep="")
      }
      contrib<-substr(contrib,1,(nchar(contrib)-2))
      nomcolum<-paste("cr",as.character(k),sep="")
      ddd[i,nomcolum]<-contrib
    }
    #add row with meaning for codes of column variables
    meaning<-""
    for (h in 1:length(cr[[k]])){
      meaning<-paste(meaning, cr[[k]][h],"('",stats[1,cr[[k]][h]],"')")
    }
    ddd[length(superan)+1,nomcolum]<-meaning
  }
  
  return(ddd)
}



## -----lasso-----
compare_lasso_agnes<-function(superan,stats) {

  k<-NULL
  for (i in 1:length(superan))
    k<-c(k,paste("C_h",i,sep=""))
  
  ddd<-data.frame(k)
  
  
  ## Graph_lasso
  library(igraph)

  ig<-as.igraph(Graph_lasso, attributes = TRUE)
  cl<-clusters(ig)
  pos<-which(cl$csize>1)

  clus<-list()
  cr<-list()

  for (i in 1:length(pos)) {
    clus[[i]]<-which(cl$membership==pos[i])
    cr[[i]]<-V(ig)$label[clus[[i]]]
  }
  
  
  for (k in 1:length(cr)) {
    for (i in 1:length(superan)) {
      print(paste("cluster",i,sep=" "))
      contrib<-""
      for (j in 1:length(cr[[k]])) {
        contrib<-paste(contrib,cr[[k]][j],"(",round(superan[[i]][cr[[k]][j]],1),"), ",sep="")
      }
      contrib<-substr(contrib,1,(nchar(contrib)-2))
      nomcolum<-paste("cr",as.character(k),sep="")
      ddd[i,nomcolum]<-contrib
    }
    #add row with meaning for codes of column variables
    meaning<-""
    for (h in 1:length(cr[[k]])){
      meaning<-paste(meaning, cr[[k]][h],"('",stats[1,cr[[k]][h]],"')")
    }
    ddd[length(superan)+1,nomcolum]<-meaning
  }
  
  return(ddd)
}


## -----Ising2-----
compare_Ising2_agnes<-function(superan,stats) {

    k<-NULL
  for (i in 1:length(superan))
    k<-c(k,paste("C_h",i,sep=""))
  
  ddd<-data.frame(k)
  

  ## Graph_Ising2
  library(igraph)

  ig<-as.igraph(Graph_Ising2, attributes = TRUE)
  cl<-clusters(ig)
  pos<-which(cl$csize>1)

  clus<-list()
  cr<-list()

  for (i in 1:length(pos)) {
    clus[[i]]<-which(cl$membership==pos[i])
    cr[[i]]<-V(ig)$label[clus[[i]]]
  }
  
  
  for (k in 1:length(cr)) {
    for (i in 1:length(superan)) {
      print(paste("cluster",i,sep=" "))
      contrib<-""
      for (j in 1:length(cr[[k]])) {
        contrib<-paste(contrib,cr[[k]][j],"(",round(superan[[i]][cr[[k]][j]],1),"), ",sep="")
      }
      contrib<-substr(contrib,1,(nchar(contrib)-2))
      nomcolum<-paste("cr",as.character(k),sep="")
      ddd[i,nomcolum]<-contrib
    }
    #add row with meaning for codes of column variables
    meaning<-""
    for (h in 1:length(cr[[k]])){
      meaning<-paste(meaning, cr[[k]][h],"('",stats[1,cr[[k]][h]],"')")
    }
    ddd[length(superan)+1,nomcolum]<-meaning
  }
  
  return(ddd)
}
