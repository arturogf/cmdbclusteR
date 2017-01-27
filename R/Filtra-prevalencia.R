# This script create new rows in the input data frame with those diagnoses 
# present in more than a threshold prevalence for our population

# -------------- define file input --------------
f = "/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-transformado-phewas.csv"
mydata = read.csv(f, header=TRUE, sep=";", fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

# -------------- define path to PHEWAS mapping file --------------
fphewas = "../data/mappings/PheWAS_code_translation_v1_2-ORIGINAL.txt"

# -------------- define file output for descriptive prevalence phewas information --------------
fphewasprev = "/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-phewas-prevalentes2.csv"

# -------------- define file output with prevalence calculations --------------
fout = "/Users/arturogf/Documents/Unidad\ Innovacion/hyponatremia/pacientes\ -\ 2011\ -\ 2015/siadh-phewas-prevalencia2.csv"

# all columns (11) are defined of type character
colClasses = c(rep("character",11))
mapeo = read.csv(fphewas, header=TRUE, sep="\t", colClasses=colClasses,fileEncoding="UTF-8", as.is = TRUE, check.names = FALSE) 

# any filtering
#mydata <- mydata[ which(mydata$Edad<=50), ]
#mydata <- mydata[ which(mydata$Edad>=18), ]

# -------------- Define the prevalence threshold (percentage) --------------
threshold_prev <- 2

# -------------- Define if the input data has been already mapped to PHEWAS --------------
input_type <- "PHEWAS" # In this case, no need to do mapping again
num_patients <- nrow(mydata)

# we suppose that first column with ICD9 separated codes is the following to Dx.Todos
pos_first_icd9 <- match("Dx.Todos", names(mydata))
pos_first_icd9 <- pos_first_icd9 + 1

# We add new row for total number of ocurrences (Total)
mydata["Total", pos_first_icd9:ncol(mydata)] <-
  colSums(mydata[, pos_first_icd9:ncol(mydata)])
# We add new row for prevalence percentage (Porcentaje)
mydata["Porcentaje", pos_first_icd9:ncol(mydata)] <-
  mydata["Total", pos_first_icd9:ncol(mydata)] / num_patients * 100
# We add a new row indicating if is prevalent (Prevalente)
mydata["Prevalente", pos_first_icd9:ncol(mydata)] <-
  (mydata["Porcentaje", pos_first_icd9:ncol(mydata)] >= threshold_prev)

# we create a filter vector for the columns (numbers) that are prevalent
filter <- integer()
for (i in pos_first_icd9:ncol(mydata))
  if (mydata["Prevalente", i] == 1)
    filter <- c(filter, i)

# those prevalent columns are merged with the patient/episode data into a new data frame 'prevalentes'
prevalentes <- mydata[, filter]
prevalentes <- cbind(mydata[, 1:pos_first_icd9 - 1], prevalentes)

# vector with the prevalent code names
codigos <- colnames(prevalentes[, pos_first_icd9:ncol(prevalentes)])

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
    phewas <- c(phewas, b)
    phewas_desc <- c(phewas_desc, c)
  }
  tabla <- data.frame(codigos, desc, phewas, phewas_desc, percent)
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
 
   # We create a file with info of phewas codes and prevalence percentages, for example:
  #"codigos";"phewas_desc";"percent"
  #"010";"Tuberculosis";3,18
  #"041.2";"Streptococcus infection";2,42 ...
  write.table(tabla,fphewasprev,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")
}

# We write a file with all the prevalence calculations (total, percentage, prevalence)
write.table(prevalentes,fout,FALSE,sep=";",row.names = FALSE, fileEncoding = "UTF-8",dec=",")


