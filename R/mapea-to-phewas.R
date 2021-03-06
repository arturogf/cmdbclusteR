# This script creates a new file with PHEWAS-coded columns as mapped from the ICD9-coded original file

# create a vector with the icd9 code strings and an empty vector for phewas codes
icd9_codes <- sort(colnames(mydata[, pos_first_ICD9:ncol(mydata)]))
phewas_codes <- vector()

# for each icd9 code we look up the phewas mapping table (column 3 is the phewas code) 
for (i in icd9_codes) {
  pos <- which(mapeo$icd9 == i)
  #--look for better method to check the following condition--
  if (length(pos) != 0) {
    newcode <- mapeo[pos, 3]
    # If not present already, we add the phewas code to the phewas_codes vector
    if (!(newcode %in% phewas_codes))
      phewas_codes <- c(phewas_codes, newcode)
  }
}

# create a new dataframe for the phewas columns and set it to 0 for the whole population
colClasses=c(rep("numeric",length(phewas_codes)))
dataph <- read.table(text = "", col.names = phewas_codes, colClasses=colClasses, as.is = TRUE, check.names = FALSE)
dataph[1:nrow(mydata),] <- 0

nomapeo<-vector()

for (l in 1:nrow(mydata)) {
  #print(l)
  # we split the codes field (string) using character ']'.
  dx_todos_sep <- unlist(strsplit(as.character(mydata[l, pos_ICD9]), '\\]'))
  for (d in dx_todos_sep) {
    #  we substring avoiding the initial '['
    h <- substring(d, 2)
    pos <- which(mapeo$icd9 == h)
    # we fill the column for the phewas code as mapped from the icd9 at pos
    if (length(pos) != 0) 
      dataph[l, mapeo[pos, 3]] <- 1
    else
      if (!(h %in% nomapeo))
        nomapeo <- c(nomapeo, h)
  }
}

# the icd9 codes without PheWAS mapping are printed out
for (i in sort(nomapeo))
  print(paste("INFO: There was no PheWAS mapping for ICD-9 code", i, sep =" "))

# Define nomapeo, where we save ICD9 codes that were not mapped to PheWAS 
nomapeo<-mydata[,sort(nomapeo)]
# and calculate how many times they are found in our data
nomapeo["total",] <- colSums(nomapeo)
nomapeo <- nomapeo["total",]
#--- code below to recover prevalent codes that were not mapped---
nomapeo["porcentaje",] <- ((nomapeo["total",])*100)/(nrow(mydata))
nomapeo <- nomapeo[(nrow(nomapeo)-1):nrow(nomapeo),]
print(paste("There was no PheWAS mapping for", ncol(nomapeo), "ICD9 codes."))

# reconsider ICD9 codes that were not mapped to PheWAS but have higher percentaje than selected prevalence threshold
remapped <- vector()
names <- vector()
for (i in 1:ncol(nomapeo)){
  if(nomapeo["porcentaje",i] >= threshold_prev){
    dataph[[names(nomapeo)[i]]] <- mydata[,names(nomapeo)[i]]
    remapped <- c(remapped, i)
    names <- c(names, names(nomapeo)[i])
    print(paste("ICD9 code", names(nomapeo)[i], 
                "has been reconsidered given its prevalence is higher than prev_threshold selected:", 
                round(nomapeo["porcentaje",i],2)))
  }
}
names(remapped) <- names
nomapeo <- nomapeo[,-remapped]
print(paste("Recovered the", length(remapped), "ICD9 codes specified above."))
#--- code above to recover prevalent codes that were not mapped---

# we can substitute ICD9 columns by PheWAS columns using mycopy defined in declaration.R
myphewas <- cbind(mycopy,dataph)

# redefine nombre_generado

print("FYI: Processed files will be saved in cmbdclusteR/data/processed by default")

nombre_generado <- paste(substring(nombre_generado, 1, nchar(nombre_generado)-4),"-phewas.csv", sep="")
fout <- file.path(directory, "data/processed", nombre_generado)
# we write the output of phewas mapping to a defined file
write.table(myphewas,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")

fout <- file.path(directory, "data/processed", paste(substring(nombre_generado,1,nchar(nombre_generado)-4),"-nomapeo.csv", sep=""))
# we write nomapeo to a defined file
write.table(nomapeo,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")

#-- remove variables that have been used, possibly using ls() (look ?remove) --
remove(dataph)

source(file.path(directory, "/R/Filtra-prevalencia.R"))