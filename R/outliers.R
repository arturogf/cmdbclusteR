detectaOutliers<-function(mij) {

#esto parece que no funciona bien
#library(DMwR)
#o1 <- outliers.ranking(dist(mij), method="sizeDiff",method.pars = NULL, clus=list(alg='hclust',meth=metodo))

library(HighDimOut)
df.fbod<-Func.FBOD(mij,10,5)
o2<-order(df.fbod,decreasing=TRUE)[1:10]

df.sod<-Func.SOD(mij,10,5,alpha = 0.8)
o3<-order(df.sod,decreasing=TRUE)[1:10]

quitar<-unique(c(o2,o3))

outliers<-mij[quitar,]
sin_outliers<-mij[-quitar,]

return(sin_outliers)

}


