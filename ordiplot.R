require(vegan)

ClusterOrdiPlot<-function(distance_measure,cluster_solution,threshold) {


groups<-names(which((table(cluster_solution))>(threshold*length(cluster_solution))))

n<-length(groups)
cols<-rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)

print(cols)
mis_etiquetas<-vector()

for (i in 1:length(groups)) {
  lab<-paste("Cluster", groups[i], sep=" ")
  mis_etiquetas<-c(mis_etiquetas,lab)
}

cmd <- cmdscale(distance_measure)
ordiplot(cmd)

#hacemos el subset de los que estan en los grupos
a<-cluster_solution %in% groups
is_in<-cluster_solution[which(a==TRUE)]

for(i in seq_along(groups)) {
  points(cmd[factor(is_in) == groups[i], ], col = cols[i], pch = 16)
}


#ordispider(cmd, factor(is_in),label = TRUE)

ordihull(cmd, factor(is_in),lty = "dotted")

legend("bottomleft", legend = mis_etiquetas, 
        pch = 16, col= cols,  bty = "n")

return(ordiplot)

}
