#Ejemplo basico con clusterboot
#convertimos el fichero con 0 y 1 a valores lógicos
data5logic<- sapply(as.data.frame(data5), as.logical)

#lo convertimos a data frame
logic5<-data.frame(data5logic)

# utilizamos el metodo clusterboot para estimar los índices de Jaccard para B ejecuciones con permutation tests
# en este caso boot (non-standard bootstrapping)

#dos opciones: k-means o hclust
km5boot <- clusterboot(logic5, B=20, bootmethod = "boot", clustermethod = kmeansCBI, krange=3,seed = 15555) 
#hemos ejecutado esta
km5boot <- clusterboot(logic5, B=100, bootmethod = "boot", clustermethod = hclustCBI, method="ward.D2", k=3) 

#creamos una matriz de distancia euclídea y hacemos el cluster jerárquico con ward.D2
d<-dist(logic5,method = "euclidean")
hcl5<-hclust(d,method="ward.D2")

#pintamos dendograma y cortamos el árbol para 3 clusters que se pintan en rojo
plot(hcl5)
groups<- cutree(hcl5, k=3)
rect.hclust(hcl5, k=3, border="red")

#La salida de los clusters es un vector con el número de cluster asignado a cada paciente
out<-capture.output(groups)
cat("hclust3", out, file="/Users/arturogf/hclust3.txt", sep=",", append=TRUE)

# Hay que encontrar y reemplazar del fichero creado con regex. \[[[:alnum:]]+\] (con espacio)
# ya que crea una especie de índices intermedios cada X pacientes. Luego asegurarse de quitar dobles espacios
# y separar por comas

