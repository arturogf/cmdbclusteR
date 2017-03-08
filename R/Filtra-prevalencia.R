# This script create new rows in the input data frame with those diagnoses 
# present in more than a threshold prevalence for our population

if(input_type == "PHEWAS"){
  mydata <- myphewas
  remove(myphewas)
}

# -------------- Define the prevalence threshold (percentage) --------------
print("Select prevalence threshold to filter you data or press ENTER for default value (2):")
threshold_prev <- as.numeric(readLines(n=1, ok=FALSE))
if (is.na(threshold_prev)){
  threshold_prev <- 2
}

num_episodes <- nrow(mydata)

# we suppose that first column with phewas separated codes is the same as the column hosting first ICD9 code
pos_first_phewas <- pos_first_ICD9

# We add new row for total number of ocurrences (Total)
mydata["Total", pos_first_phewas:ncol(mydata)] <-
  colSums(mydata[, pos_first_phewas:ncol(mydata)])
# We add new row for prevalence percentage (Porcentaje)
mydata["Porcentaje", pos_first_phewas:ncol(mydata)] <-
  mydata["Total", pos_first_phewas:ncol(mydata)] / num_episodes * 100
# We add a new row indicating if is prevalent (Prevalente)
mydata["Prevalente", pos_first_phewas:ncol(mydata)] <-
  (mydata["Porcentaje", pos_first_phewas:ncol(mydata)] >= threshold_prev)

# we create a filter vector for the columns (numbers) that are prevalent
filter <- integer()
for (i in pos_first_phewas:ncol(mydata))
  if (mydata["Prevalente", i] == 1)
    filter <- c(filter, i)

#-- check if we can avoid creating a new dataframe 'prevalentes' --

# those prevalent columns are merged with the patient/episode data into a new data frame 'prevalentes'
prevalentes <- mydata[, filter]
prevalentes <- cbind(mydata[, 1:pos_first_phewas - 1], prevalentes)

# vector with the prevalent code names
codigos <- colnames(prevalentes[, pos_first_phewas:ncol(prevalentes)])

# these vectors are used to create a descriptive output table of prevalent phewas-encoded diagnoses
desc <- vector()
phewas <- vector()
phewas_desc <- vector()
percent <- vector()

# If file is not PHEWAS
if (!(input_type == "PHEWAS")) {
  for (i in codigos) {
    percent <- c(percent, prevalentes["Porcentaje", i])
    pos <- which(mapeo$icd9 == i)
    if (length(pos) == 0) {
      a <- "NO DISPONIBLE"
      b <- "NO DISPONIBLE"
      c <- "NO DISPONIBLE"
    }
    else {
      #select into a,b,c the columns number 2 (description), 3 (phewas), and 4 (desc-phewas) of mapping file
      a <- mapeo[pos, 2]
      b <- mapeo[pos, 3]
      c <- mapeo[pos, 4]
    }
    desc <- c(desc, a)
    #phewas <- c(phewas, b)
    #phewas_desc <- c(phewas_desc, c)
  }
  #tabla <- data.frame(codigos, desc, phewas, phewas_desc, percent)
  tabla <- data.frame(codigos, desc, percent)
  
    # pos <- 1
  # tabla[["phewas_percent"]] <- rep(0,nrow(tabla))
  # for(i in 1:nrow(tabla)){
  #   if(tabla[i,"phewas"] == tabla[pos,"phewas"]){
  #     tabla[pos, "phewas_percent"] <-tabla[pos, "phewas_percent"] + tabla[i,"percent"] 
  #   }else{
  #     pos <- i
  #     tabla[pos, "phewas_percent"] <-tabla[pos, "phewas_percent"] + tabla[i,"percent"]
  #   }
  # }
}

# If file is already PHEWAS
if (input_type == "PHEWAS") {
  for (i in codigos) {
    percent <- c(percent, prevalentes["Porcentaje", i])
    pos <- which(mapeo$phewas_code == i)
    if (length(pos) == 0)
      c <- "NO DISPONIBLE"
    else
      c <- mapeo[pos[1], 4] #select fourth column (desc-phewas) from mapping file
    phewas_desc <- c(phewas_desc, c)
  }
  tabla <- data.frame(codigos, phewas_desc, percent)
  tabla[order(tabla[["percent"]], decreasing=TRUE),]
}

print("FYI: Processed file will be saved in cmbdclusteR/data/processed by default")
# ------ define file output for descriptive prevalence phewas information --------------
# redefine nombre_generado
nombre_generado <- paste(substring(nombre_generado, 1, nchar(nombre_generado)-4),"-prevalence", threshold_prev, ".csv", sep="")
fphewasprev <- file.path(directory, "data/processed", nombre_generado)

# We create a file with info of phewas codes and prevalence percentages, for example:
#"codigos";"phewas_desc";"percent"
#"010";"Tuberculosis";3,18
#"041.2";"Streptococcus infection";2,42 ...
write.table(tabla,fphewasprev,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")

# ------ define file output with prevalence calculations --------------
#redefine nombre_generado
nombre_generado <- paste(substring(nombre_generado, 1, nchar(nombre_generado)-4),"Filtered.csv", sep="")
fout <- file.path(directory, "data/processed", nombre_generado)

# We write a file with all the prevalence calculations (total, percentage, prevalence)
write.table(prevalentes,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")

remove(tabla)
source(file.path(directory, "/R/asignaclusters.R"))