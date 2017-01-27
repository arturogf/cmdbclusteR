# Este programa crea nuevas filas en el data frame de entrada
# con aquellas patologías presentes en más del 5 por ciento
f="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-transformado-phewas.csv"
mydata=read.csv(f, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

#mydata <- mydata[ which(mydata$Edad<=50), ]
#mydata <- mydata[ which(mydata$Edad>=18), ]

# Definir nivel de prevalencia deseado
threshold_prev<-2
tipo_entrada<-"PHEWAS" # significa que lo que entra ya tiene codigos phewas y no hay que hacer mapeo
num_patients<-nrow(mydata)

# suponemos que la primera columna con codigos es la siguiente a Dx.Todos
pos_primer_icd9<-match("Dx.Todos",names(mydata))
pos_primer_icd9<-pos_primer_icd9 + 1

# Añadimos las filas de total, porcentaje y 
mydata["Total",pos_primer_icd9:ncol(mydata)]<-colSums(mydata[,pos_primer_icd9:ncol(mydata)])
mydata["Porcentaje",pos_primer_icd9:ncol(mydata)]<-mydata["Total",pos_primer_icd9:ncol(mydata)]/num_patients*100
mydata["Prevalente",pos_primer_icd9:ncol(mydata)]<-(mydata["Porcentaje",pos_primer_icd9:ncol(mydata)]>=threshold_prev)

# creamos un vector con las posiciones que tienen la última fila "Prevalente" a 1
filter<-integer()
for (i in pos_primer_icd9:ncol(mydata)) 
	if (mydata["Prevalente",i]==1) 
		filter<-c(filter,i)

# seleccionamos esas columnas y las juntamos con las primeras (13) columnas
prevalentes<-mydata[,filter]
prevalentes<-cbind(mydata[,1:pos_primer_icd9-1],prevalentes)

# Defina el path al fichero de mapeo phewas
fphewas="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/Henar/PheWAS_code_translation_v1_2-ORIGINAL.txt"
# todas las columnas (11) son de tipo character
colClasses=c(rep("character",11))
mapeo=read.csv(fphewas, header=TRUE, sep="\t", colClasses=colClasses,fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

#definimos vectores que usaremos para luego ponerlos de columnas en el data frame de prevalencias
desc<-vector()
phewas<-vector()
phewas_desc<-vector()
percent<-vector()

codigos<-colnames(prevalentes[,pos_primer_icd9:ncol(prevalentes)])

if (!(tipo_entrada=="PHEWAS")) {
	for (i in codigos) { 
		percent<-c(percent,prevalentes["Porcentaje",i])
		pos<-which(mapeo$icd9==i)
		if (length(pos)==0) {
			a<-"NO DISPONIBLE"
			b<-"NO DISPONIBLE"
			c<-"NO DISPONIBLE"
		}
		else { #seleccionamos la columna 2 (descripcion) y 3 (phewas), 4 (desc-phewas) del fichero de mapeo
			a<-mapeo[pos,2]	
			b<-mapeo[pos,3]
			c<-mapeo[pos,4]
		}
		desc<-c(desc,a)
		phewas<-c(phewas,b)
		phewas_desc<-c(phewas_desc,c)
	}
	tabla<-data.frame(codigos,desc,phewas,phewas_desc,percent)
}
if (tipo_entrada=="PHEWAS") { 
	for (i in codigos) { 
		percent<-c(percent,prevalentes["Porcentaje",i])
		pos<-which(mapeo$phewas_code==i)
		if (length(pos)==0) 
			c<-"NO DISPONIBLE"
		else
			c<-mapeo[pos[1],4] #seleccionamos la columna 4 (desc-phewas) del fichero de mapeo
		phewas_desc<-c(phewas_desc,c)
	}
	tabla<-data.frame(codigos,phewas_desc,percent)
}

# salida a fichero salida de la tabla de estadisticas de prevalentes
fout="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-phewas-prevalentes2.csv"
write.table(tabla,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")

# salida a fichero salida de los calculos completos de prevalencia
fout="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-phewas-prevalencia2.csv"
write.table(prevalentes,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")


