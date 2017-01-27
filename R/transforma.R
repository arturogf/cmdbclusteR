# Este programa recibe un fichero csv (por ;) del CMBD con los siguientes campos:
# Año	
# Tipo Actividad	
# Nº Historia	
# Sexo -Descripción-	
# Fecha Nacimiento	
# Edad	
# Fecha Ingreso	
# Fecha Alta	
# Sección Alta
# Motivo Alta -Descripción-	
# GRD	
# Versión	
# Dx Todos
# 
# Escribe como salida un fichero "feature matrix" donde cada código de Dx.Todos
# se añade a una columna y se rellena con 0/1 para cada paciente.

# Definimos el numero de columna donde estan los dx
coldx<-match("Dx.Todos",names(mydata))

# definir el fichero csv a tratar. No olvidar que esté en UTF-8 y con UNIX (LF)

#f<- file.choose()
f="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/Hiponatremia-SIADH-FechaNac.csv"
fout="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-mayor18-transformado.csv"
mydata=read.csv(f, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE) 

# los mayores de 18
mydata <- mydata[ which(mydata$Edad>18), ]

# si queremos por ejemplo los que tienen siadh
mydata<-mydata[grep("253.6",mydata[[coldx]],invert = FALSE),]

# Copiamos todos los valores de celdas superiores para las celdas vacías
# esto se debe a la agrupación por fecha que hacen los de codificación en el excel,
# al exportar a csv se pierden valores de las celdas agrupadas.
for (i in 1:nrow(mydata)) 
	for (j in 1:match("Fecha.Alta",names(mydata))) 
		if (is.na(mydata[i,j]) | mydata[i,j]=="") {
			print(paste(as.character(i),as.character(j),sep=","))
			if (i>1)
				mydata[i,j]<-mydata[i-1,j]
		}

#hacemos una copia de los campos iniciales 
mycopy<-mydata

# añadimos las columnas de los ICD-9 existentes en el conjunto de datos
for (l in 1:nrow(mydata)) {
	# hacemos un split de Dx.Todos por el caracter ']'
	dx_todos_sep<-unlist(strsplit(as.character(mydata[l,coldx]),'\\]'))
	for (d in dx_todos_sep) {
		# quitamos el corchete de inicio 
		h<-substring(d,2)
		if(h %in% colnames(mydata))
		{
			# ponemos un 1 para ese código en esa fila
		 	mydata[l,match(h,names(mydata))]<-1
		}
		else
			{
				 #añadir una columna con todo a 0 para ese ICD-9
				 print(paste("Añadiendo código ICD-9 ", h))
				 mydata[[h]]<-0	 
				 # ponemos un 1 para ese código en esa fila
				 mydata[l,match(h,names(mydata))]<-1
			}
	}
}

# hacemos el subset de los nuevos campos añadidos y los ordenamos
pos_primer_icd9<-coldx + 1
mydata<-mydata[c(pos_primer_icd9:ncol(mydata))]
mydata<-mydata[ , order(names(mydata))]

# unimos los dos data frames en uno
mezcla<-cbind(mycopy,mydata)

# salida a fichero salida
write.table(mezcla,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8")