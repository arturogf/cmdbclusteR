# The scripts outputs a feature matrix file where each 
# of the ICD9 codes is added as a column and is filled with 0/1 
# for each patient.

# For each unique ICD9 existing code, we add a new column. Processed by row.
for (l in 1:nrow(mydata)) {
  # we split the ICD9 codes field (string) using character ']'.
  dx_todos_sep <- unlist(strsplit(as.character(mydata[l, pos_ICD9]), '\\]'))
  for (d in dx_todos_sep) {
    # we substring avoiding the initial '['
    h <- substring(d, 2)
    # if the ICD9 code has been already added as column, we indicate the occurrence. *improve code*
    if (h %in% colnames(mydata))
    {
      mydata[l, match(h, colnames(mydata))] <- 1
    }
    # Otherwise, add column and indicate the ocurrence.
    else
    {
      print(paste("INFO: Adding column for ICD-9 code ", h))
      mydata[[h]] <- 0
      mydata[l, match(h, colnames(mydata))] <- 1
    }
  }
}

# Reorder ICDs columns by name.
icd9_codes <- sort(colnames(mydata[,c(pos_first_ICD9:ncol(mydata))]))
mydata <- mydata[,icd9_codes]
mydata <- cbind(mycopy, mydata) 

# -------------- define file output --------------
nombre_generado <-
  paste(
    tools::file_path_sans_ext(basename(f)),
    diagnosefilter,
    GRDfilter,
    agefilter,
    genderfilter,
    "transformado.csv",
    sep = "-"
  )

#create a file to save processed data
if (!file.exists(file.path(directory, "data/processed"))){
  dir.create(file.path(directory, "data/processed"))
}

print("FYI: Processed file will be saved in cmbdclusteR/data/processed by default")
fout <- file.path(directory, "data/processed", nombre_generado)

# write the data frame into the output file, a ; separated CSV that can be directly imported from excel
write.table(mydata, fout, FALSE, sep=";", row.names = FALSE)

# Select whether you want to map data to PHEWAS or not
print("Write Y to map data to PheWAS code or press ENTER to continue:")
optionmap <- readLines(n=1, ok=FALSE)

if (optionmap == "Y"){
  input_type <- "PHEWAS"
  source(file.path(directory, "/R/mapea-to-phewas.R"))
}else{
  input_type <- "GIVEN"
  source(file.path(directory, "/R/Filtra-prevalencia.R"))
}