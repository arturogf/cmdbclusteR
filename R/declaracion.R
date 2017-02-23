# -------------- define path to file input --------------
print("A window will pop up to select data input file.")
f <- tk_choose.files(default = "", caption = "Select data input file", multi = FALSE)
# Read data separated by ";"
mydata = read.csv(f, header=TRUE, sep=";", as.is = TRUE) 

# -------------- define path to mapping file --------------
print("A window will pop up to select mapping file (ICD9-PheWAS or ICD10-PheWAS).")
fmapping <- tk_choose.files(default = "", caption = "Select mapping file", multi = FALSE)
# Read data separated by "\t"
mapping = read.csv(fmapping, header=TRUE, sep="\t", as.is = TRUE, check.names=FALSE)
for (i in 1:ncol(mapping)){
  class(mapping[[i]]) <- "character"
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

# Define the number of column that hosts patients' age
print("Number of column that hosts patients' ages (write 0 if there is no such column):")
pos_age <- as.numeric(readLines(n=1, ok=FALSE))
if(pos_age==0){
  print("Number of column that hosts patients' dates of birth:")
  pos_birthdate <- as.numeric(readLines(n=1, ok=FALSE))
  mydata[,pos_birthdate] <- gsub("/", "-", mydata[,pos_birthdate])
  # Create new column with patients' ages
  mydata$Edad <- round(as.numeric(as.Date(mydata[[pos_discharge]], "%d-%m-%Y")-as.Date(mydata[[pos_birthdate]], "%d-%m-%Y"))/365)
  pos_age <- ncol(mydata)
}
# Define the number of column that hosts patients' gender
print("Number of column that hosts patients' gender:")
pos_gender <- as.numeric(readLines(n=1, ok=FALSE))

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
  print("Choose up threshold:")
  umbral2 <- as.numeric(readLines(n=1, ok=FALSE))
  mydata <- mydata[ which(mydata[[pos_age]]>umbral1 & mydata[[pos_age]]<umbral2), ]
  agefilter <- paste("between", umbral1, "and", umbral2, sep="")
}else if(intervalo==""){
  agefilter<-"allages"
}

#prueba
#Or choose a gender filter
print("Choose a gender filter: write M for MALE, F for FEMALE and B for BOTH:")
gender <- readLines(n=1, ok=FALSE)
if(gender == "M"){
  mydata <- mydata[ which(mydata[[pos_gender]]=="Hombre" | mydata[[pos_gender]]==1), ]
  genderfilter <- "male"
}else if(gender == "F"){
  mydata <- mydata[ which(mydata[[pos_gender]]=="Mujer" | mydata[[pos_gender]]==2), ]
  genderfilter <- "female"
}else if(gender == "B" | gender==""){
  genderfilter <- "both"
}

# Select GRD filter
if(pos_GRD!=0){
  print("Select GRD filter (if filtering by several GRD values, separate them by ','): ")
  GRDvalues <- as.numeric(unlist(strsplit(readLines(n=1,ok=FALSE), split=",")))
  GRDchar <- paste(as.character(sort(GRDvalues)), collapse="-")
  GRDfilter <- paste("GRD", GRDchar, sep="-")
  mydata <- mydata[which(mydata[[pos_GRD]] %in% GRDvalues),]
}

#SIADH ICD9: 253.6
#mydata <- mydata[grep("253.6", mydata[[pos_ICD9]], invert = FALSE),]
#diagnosefilter <- "siadh"

# or only does that present some diagnose, e. g. humour disorder
mydata <- mydata[grep("296", mydata[[pos_ICD9]], invert = FALSE),]
diagnosefilter <- "humour"

pos_first_ICD9 <- ncol(mydata) + 1