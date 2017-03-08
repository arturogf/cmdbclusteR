answer <- ""
if (file.exists(file.path(directory, "setup.cfg"))){
  load(file.path(directory, "setup.cfg"))
  print(paste("You are using", basename(f), "input data. Write F to select another input data or press ENTER to continue:"))
  answer <- as.character(readLines(n=1, ok=FALSE))
  if(answer==""){
    mydata = read.csv(f, header=TRUE, sep=";", as.is = TRUE)
    mydata[,pos_entry] <- gsub("/", "-", mydata[,pos_entry])
    mydata[,pos_discharge] <- gsub("/", "-", mydata[,pos_discharge])
    mapeo = read.csv(fmapping, header=TRUE, sep="\t", encoding="UTF-8", as.is = TRUE, check.names=FALSE)
    for (i in 1:ncol(mapeo)){
      class(mapeo[[i]]) <- "character"
    }
  }
}
if(!(file.exists(file.path(directory, "setup.cfg"))) | answer=="F"){
  # -------------- define path to file input --------------
  print("A window will pop up to select data input file.")
  f <- tk_choose.files(default = "", caption = "Select data input file", multi = FALSE)
  # Read data separated by ";"
  mydata = read.csv(f, header=TRUE, sep=";", as.is = TRUE) 
  
  # -------------- define path to mapping file --------------
  print("A window will pop up to select mapping file (ICD9-PheWAS or ICD10-PheWAS).")
  fmapping <- tk_choose.files(default = "", caption = "Select mapping file", multi = FALSE)
  # Read data separated by "\t"
  mapeo = read.csv(fmapping, header=TRUE, sep="\t", as.is = TRUE, check.names=FALSE)
  for (i in 1:ncol(mapeo)){
    class(mapeo[[i]]) <- "character"
  }
  
  # Define the number of column that hosts the patient ID / NHC
  print("Number of column that hosts the NHC:")
  pos_numeroHC <- as.numeric(readLines(n=1, ok=FALSE))
  
  # Define the number of column that hosts the ICD9 codes
  print("Number of column that hosts the ICD9 codes:")
  pos_ICD9 <- as.numeric(readLines(n=1, ok=FALSE))
  
  # Define the number of column that hosts the GRD
  print("Number of column that hosts GRD (write 0 if there is no such column):")
  pos_GRD <- as.numeric(readLines(n=1, ok=FALSE))
  
  # Define the number of column that hosts the entry date
  print("Number of column that hosts patients' entry date:")
  pos_entry <- as.numeric(readLines(n=1, ok=FALSE))
  mydata[,pos_entry] <- gsub("/", "-", mydata[,pos_entry])
  
  # Define the number of column that hosts the discharge date
  print("Number of column that hosts patients' discharge date:")
  pos_discharge <- as.numeric(readLines(n=1, ok=FALSE))
  mydata[,pos_discharge] <- gsub("/", "-", mydata[,pos_discharge])
  
  # Define the number of column that hosts patients' gender
  print("Number of column that hosts patients' gender:")
  pos_gender <- as.numeric(readLines(n=1, ok=FALSE))
  
  # Define the number of column that hosts patients' age
  print("Number of column that hosts patients' ages (write 0 if there is no such column):")
  pos_age <- as.numeric(readLines(n=1, ok=FALSE))
  
  if(pos_age==0){
    print("Number of column that hosts patients' dates of birth:")
    pos_birthdate <- as.numeric(readLines(n=1, ok=FALSE))
  }else{
    pos_birthdate <- integer()
  }
  
  # Save column positions and file names for this data input
  save(f, fmapping, pos_age, pos_birthdate, pos_discharge, pos_entry, pos_gender, pos_GRD, pos_ICD9, pos_numeroHC, file=file.path(directory, "setup.cfg"))
}

# Create a new column age if pos_age equals 0
if(pos_age==0){
  mydata[,pos_birthdate] <- gsub("/", "-", mydata[,pos_birthdate])
  # Create new column with patients' ages
  mydata$Edad <- round(as.numeric(as.Date(mydata[[pos_discharge]], "%d-%m-%Y")-as.Date(mydata[[pos_birthdate]], "%d-%m-%Y"))/365)
  pos_age <- ncol(mydata)
}


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

# Create new column with patient's stay in days
mydata$Estancia <- round(as.numeric(as.Date(mydata[[pos_discharge]], "%d-%m-%Y")-as.Date(mydata[[pos_entry]], "%d-%m-%Y")))
pos_stay <- ncol(mydata)

# Add any kind of filtering that we want on the data, e.g. age > 18
print("Choose an age filter (or press ENTER to select all ages): write L for LESS THAN, H for HIGHER THAN or B for BETWEEN:")
intervalo <- readLines(n=1, ok=FALSE)
if(intervalo == "L"){
 print("Choose top threshold:")
 umbral <- as.numeric(readLines(n=1, ok=FALSE))
 mydata <- mydata[ which(mydata[[pos_age]]<umbral), ]
 agefilter <- paste("lessthan", umbral, sep="")
}else if(intervalo == "H"){
 print("Choose down threshold:")
 umbral <- as.numeric(readLines(n=1, ok=FALSE))
 mydata <- mydata[ which(mydata[[pos_age]]>umbral), ]
 agefilter <- paste("higherthan", umbral, sep="")
}else if(intervalo == "B"){
 print("Choose down threshold:")
 umbral1 <- as.numeric(readLines(n=1, ok=FALSE))
 print("Choose top threshold:")
 umbral2 <- as.numeric(readLines(n=1, ok=FALSE))
 mydata <- mydata[ which(mydata[[pos_age]]>umbral1 & mydata[[pos_age]]<umbral2), ]
 agefilter <- paste("between", umbral1, "and", umbral2, sep="")
}else if(intervalo==""){
 agefilter<-"allages"
}

#Or choose a gender filter
print("Choose a gender filter: write M for MALE, F for FEMALE and B for BOTH:")
gender <- readLines(n=1, ok=FALSE)
if(gender == "M"){
 mydata <- mydata[ which(tolower(mydata[[pos_gender]])=="hombre" | mydata[[pos_gender]]==1), ]
genderfilter <- "male"
}else if(gender == "F"){
 mydata <- mydata[ which(tolower(mydata[[pos_gender]])=="mujer" | mydata[[pos_gender]]==2), ]
 genderfilter <- "female"
}else{
 genderfilter <- "both"
}

# Select GRD filter
if(pos_GRD!=0) {
  print("Select GRD filter (if filtering by several GRD values, separate them by ','): ")
  GRDvalues <-
    trimws(unlist(strsplit(readLines(n = 1, ok = FALSE), split = ",")))
  if (length(GRDvalues) == 0) {
    GRDfilter <- "allGRD"
  } else{
    GRDchar <-
      paste(as.character(sort(as.numeric(GRDvalues))), collapse = "-")
    GRDfilter <- paste("GRD", GRDchar, sep = "-")
  }
  matches <-
    unique(grep(paste(GRDvalues, collapse = "|"), mydata[[pos_GRD]]))
  mydata <- mydata[matches, ]
} else{
  GRDfilter<-""
}

#Select ICD9 filter
print("Choose 1 for filtering just ICD9 found in D1 (primary diagnose) or 2 for ICD9 found in all diagnoses:")
filtermode <- readLines(n=1, ok=FALSE)
print("Select ICD9 filter (if filtering by several ICD9 values, separate them by ','): ")
ICD9values <- trimws(unlist(strsplit(readLines(n=1,ok=FALSE), split=",")))
if(filtermode=="1"){
 if(length(ICD9values)==0){
   diagnosefilter <- "allICD9"
 }else{
   ICD9char <- paste(as.character(sort(as.character(ICD9values))), collapse="-")
   diagnosefilter <- paste("ICD9", ICD9char, sep="-")
 }
 d1 <- rep(0, nrow(mydata))
 for (i in 1:nrow(mydata)){
   d1[i] <- substring(mydata[i,pos_ICD9],2,regexpr("]", mydata[i,pos_ICD9])-1)
 }
 matchesICD9 <- unique(grep(paste(ICD9values,collapse="|"), d1))
 mydata <- mydata[matchesICD9,]
 filtermode<-"d1"
}else{
 if(length(ICD9values)==0){
   diagnosefilter <- "allICD9"
 }else{
   ICD9char <- paste(as.character(sort(as.character(ICD9values))), collapse="-")
   diagnosefilter <- paste("ICD9", ICD9char, sep="-")
 }
 matchesICD9 <- unique(grep(paste(ICD9values,collapse="|"), mydata[[pos_ICD9]]))
 mydata <- mydata[matchesICD9,]
 filtermode <- "dall"
}

mycopy <- mydata
pos_first_ICD9 <- ncol(mydata) + 1

source(file.path(directory,"/R/transforma.R"))