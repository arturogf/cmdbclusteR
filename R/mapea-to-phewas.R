# Este programa crea un nuevo fichero con codigos phewas
#desde la población de individuos originales
f="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-mayor18-transformado.csv"
mydata=read.csv(f, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

# Defina el path al fichero de mapeo phewas
fphewas="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/Henar/PheWAS_code_translation_v1_2-ORIGINAL.txt"
# todas las columnas (11) son de tipo character
colClasses=c(rep("character",11))
mapeo=read.csv(fphewas, header=TRUE, sep="\t", colClasses=colClasses,fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

# suponemos que la primera columna con codigos es la siguiente a Dx.Todos
pos_primer_icd9<-match("Dx.Todos",names(mydata))
pos_primer_icd9<-pos_primer_icd9 + 1

codigos<-colnames(mydata[,pos_primer_icd9:ncol(mydata)])

phewas_codes<-vector()

# miramos a que mapean todos los icd9
for (i in codigos) { 
	pos<-which(mapeo$icd9==i)
	if (length(pos)!=0) {
		newcode<-mapeo[pos,3]
		# añadimos al vector el código phewas
		if (!(newcode %in% phewas_codes))
			phewas_codes<-c(phewas_codes,newcode)
	}
}

colClasses=c(rep("numeric",length(phewas_codes)))

dataph<-read.table(text = "",col.names = phewas_codes, colClasses=colClasses, as.is = TRUE, check.names = FALSE)
dataph[1:nrow(mydata),]<-0


nomapeo<-vector()

for (l in 1:nrow(mydata)) {
	#print(l)
	# hacemos un split por el caracter ']'
	dx_todos_sep<-unlist(strsplit(as.character(mydata[l,match("Dx.Todos",names(mydata))]),'\\]'))
	for (d in dx_todos_sep) {
		# quitamos el corchete de inicio 
		h<-substring(d,2)
		pos<-which(mapeo$icd9==h)
		if (length(pos)!=0)
		 	dataph[l,mapeo[pos,3]]<-1
		 else 
		 	if (!(h %in% nomapeo))
		 		nomapeo<-c(nomapeo,h)
	}
}

for (i in nomapeo)
	print(paste("No se encontró mapeo del ICD-9",i,sep=" "))

# copiamos las primeras columnas
pos_primer_icd9<-match("Dx.Todos",names(mydata))
mycopy<-mydata[,c(1:pos_primer_icd9)]

# unimos los dos data frames en uno
mezcla<-cbind(mycopy,dataph)

# salida a fichero salida de la tabla de estadisticas de prevalentes
fout="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-mayor18-transformado-phewas.csv"
write.table(mezcla,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")
