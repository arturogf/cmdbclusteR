# This script creates a new file with PHEWAS-coded columns as mapped from the ICD9-coded original file

# -------------- define path to PHEWAS mapping file --------------
fphewas="/Users/arturogf/cmdbclusteR/data/mappings/PheWAS_code_translation_v1_2-ORIGINAL.txt"
# indicate that all columns (11) are of type character 
colClasses=c(rep("character",11))
mapeo=read.csv(fphewas, header=TRUE, sep="\t", colClasses=colClasses,fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

# -------------- define file input --------------
f="/Users/arturogf/cmdbclusteR/data/processed/siadh-mayor18-transformado.csv"
mydata=read.csv(f, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

# -------------- define file output --------------
fout="/Users/arturogf/cmdbclusteR/data/processed/siadh-mayor18-transformado-phewas.csv"

# we suppose that first column with ICD9 separated codes is the following to Dx.Todos
pos_first_icd9 <- match("Dx.Todos", names(mydata))
pos_first_icd9 <- pos_first_icd9 + 1

# create a vector with the icd9 code strings and an empty vector for phewas codes
icd9_codes <- colnames(mydata[, pos_first_icd9:ncol(mydata)])
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

colClasses=c(rep("numeric",length(phewas_codes)))

# create a new dataframe for the phewas columns and set it to 0 for the whole population
dataph <- read.table(text = "",col.names = phewas_codes, colClasses=colClasses, as.is = TRUE, check.names = FALSE)
dataph[1:nrow(mydata),] <- 0

nomapeo<-vector()

for (l in 1:nrow(mydata)) {
  #print(l)
  # we split the codes field (string) using character ']'.
  dx_todos_sep <-
    unlist(strsplit(as.character(mydata[l, match("Dx.Todos", names(mydata))]), '\\]'))
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
for (i in nomapeo)
  print(paste("INFO: There was no PheWAS mapping for ICD-9 code", i, sep =" "))

#--in this case we can avoid the mycopy by combining mydata+dataph and remove the icd-9 separate columns--

# we do a copy of the first columns before the icd9 codes 
pos_first_icd9 <- match("Dx.Todos", names(mydata))
mycopy <- mydata[, c(1:pos_first_icd9)]

# merge the two data frames into one
mezcla <- cbind(mycopy,dataph)

# we write the output of phewas mapping to a defined file
write.table(mezcla,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")

#-- remove variables that have been used, possibly using ls() (look ?remove) --
remove(dataph)
remove(mapeo)
