# cmdbclusteR
Some R scripts for clustering patient data from patient discharge data from Spain databases (CMBD). 

The code would need to be refactored properly with functions.  Scripts usually takes a file as input and create new files as output, with a proper workflow that takes from the original CMBD file to the output csv file with the clustering results.

The scientific workflow works as follows (ordered by sequence):

- transforma.R -> input the base CSV file with ICD9 and outputs a feature matrix file with separated ICD9-columns
- mapea-to-phewas.R -> input the ICD9 feature matrix file and outputs a PHEWAS feature matrix file
- Filtra-prevalencia.R -> it filters out the PHEWAS columns not reaching a stablished prevalence threshold for the dataset 
- asignaclusters.R -> this script receives the prevalence-filtered file and add to it a new column with the assigned cluster for each row. It uses hierarchical clustering and similarity profile analysis for the task.