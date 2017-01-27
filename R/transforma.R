# This script receives a CSV file with the following fields
# Year (Año)
# Type of Activity (Tipo de Actividad)	
# Patient Record (Num. Historia Clinica)
# Gender (Sexo)	
# Birth date (Fecha Nacimiento)
# Age (Edad)
# Admission Date (Fecha Ingreso)
# Discharge Date (Fecha Alta)	
# Discharge Service (Sección Alta)
# Discharge Reason (Motivo Alta)	
# DRG (GRD)
# Version (Versión)	
# ICD9 Diagnosis Codes (Dx Todos)
# 
# The scripts outputs a feature matrix file where each 
# of the ICD9 codes is added as a column and is filled with 0/1 
# for each patient.

# Files must be in UTF-8 encoding, UNIX (LF)
#f<- file.choose()
# -------------- define file input --------------
f="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/Hiponatremia-SIADH-FechaNac.csv"
# -------------- define file output --------------
fout="/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-mayor18-transformado.csv"

# Read data separated by ';' by default when exporting from excel
mydata=read.csv(f, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE) 

# Define the number of column that host the ICD9 codes
coldx <- match("Dx.Todos",names(mydata))

# Add any kind of filtering that we want on the data, e.g. age > 18
mydata <- mydata[ which(mydata$Edad>18), ]

# or only does that present SIADH
mydata<-mydata[grep("253.6",mydata[[coldx]],invert = FALSE),]

# --------- THIS LOOP IS ONLY NEEDED IN OUR SETTINGS ----------
# For each empty cell, we fill it with its upper cell value. This is due to the source data
# having been grouped by date. When exporting excel data to csv, the lower grouped data values are lost.
for (i in 1:nrow(mydata))
  for (j in 1:match("Fecha.Alta", names(mydata)))
    if (is.na(mydata[i, j]) | mydata[i, j] == "") {
      print(paste(as.character(i), as.character(j), sep = ","))
      if (i > 1)
        mydata[i, j] <- mydata[i-1, j]
    }

# made a copy of the initial data fields
mycopy <- mydata

# For each unique ICD9 existing code, we add a new column. Processed by row.
for (l in 1:nrow(mydata)) {
  # we split the ICD9 codes field (string) using character ']'.
  dx_todos_sep <- unlist(strsplit(as.character(mydata[l, coldx]), '\\]'))
  for (d in dx_todos_sep) {
    # we substring avoiding the initial '['
    h <- substring(d, 2)
    # if the ICD9 code has been already added as column, we indicate the occurrence. *improve code*
    if (h %in% colnames(mydata))
    {
      mydata[l, match(h, names(mydata))] <- 1
    }
    # Otherwise, add column and indicate the ocurrence.
    else
    {
      print(paste("INFO: Adding column for ICD-9 code ", h))
      mydata[[h]] <- 0
      mydata[l, match(h, names(mydata))] <- 1
    }
  }
}

# we put in mydata only the new added columns and reorder columns by name.
pos_first_icd9 <- coldx + 1
mydata <- mydata[c(pos_first_icd9:ncol(mydata))]
mydata <- mydata[, order(names(mydata))]

# we join the two data frames into one
mezcla <- cbind(mycopy, mydata)

# write the data frame into the output file, a ; separated CSV that can be directly imported from excel
write.table(mezcla,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8")