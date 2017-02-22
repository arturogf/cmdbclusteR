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

# --------- THIS LOOP IS ONLY NEEDED IN OUR SETTINGS ----------
# For each empty cell, we fill it with its upper cell value. This is due to the source data
# having been grouped by date. When exporting excel data to csv, the lower grouped data values are lost.
for (i in 1:nrow(mydata))
  for (j in 1:pos_discharge)
    if (is.na(mydata[i, j]) | mydata[i, j] == "") {
      print(paste(as.character(i), as.character(j), sep = ","))
      if (i > 1)
        mydata[i, j] <- mydata[i-1, j]
    }

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

# we reorder ICDs columns by name.##REVISAR!
first_icd9 <- pos_ICD9 + 1
ordered <- mydata[c(first_icd9:ncol(mydata))]
ordered <- mydata[, order(names(mydata))]
mydata <- cbind(mydata[c(1:(first_icd9-1))],ordered) 


# -------------- define file output --------------
nombre_generado <-
  paste(
    tools::file_path_sans_ext(basename(f)),
    diagnosefilter,
    agefilter,
    genderfilter,
    "transformado.csv",
    sep = "-"
  )

print("FYI: Processed file will be saved in cmbdclusteR-master/data/processed by default")
fout <- file.path(directory, "data/processed", nombre_generado)

# create a file to save processed data
if (!file.exists(file.path(directory, "data/processed"))){
  dir.create(file.path(directory, "data/processed"))
}

# write the data frame into the output file, a ; separated CSV that can be directly imported from excel
write.table(mydata, fout, FALSE, sep=";", row.names = FALSE)

#source(file.path(directory, "/R/mapea-to-phewas.R"))