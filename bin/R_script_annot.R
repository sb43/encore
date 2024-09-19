##This is a Rcript to annotate the merge counts file with the information comming from the library
##The first argument corresponds to the location of the all.counts file you want to annotate
##The second argument corresponds to the location of the _index.sorted.txt file with the library anotation 
##The third argument is the name and location of the output file. 
# Check if enough arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) < 3) {
  stop("Insufficient arguments. Please provide input, library and output file paths.")
}

# Get arguments from command line arguments
input_file <- commandArgs(trailingOnly = TRUE)[1]
lib_file <- commandArgs(trailingOnly = TRUE)[2]
output_file <- commandArgs(trailingOnly = TRUE)[3]

# Check if the input file  and the lib file exists
if (!file.exists(input_file)) {
  stop("Input file does not exist.")
}
if (!file.exists(lib_file)) {
  stop("Lib file does not exist.")
}

# Read the files 
data <- as.matrix(read.delim(input_file))
libs <- as.matrix(read.delim(lib_file))

###Using rownames to join, we leave the first column behind

rownames(data)=data[,'sgRNA']
libs=cbind(libs,data[libs[,'ID'],2:ncol(data)])

# Save the data to the specified output file
write.table(libs, file = output_file, row.names=F, quote = F,sep = "\t")

# Print a message indicating successful execution
cat("TXT file saved successfully at", output_file, "\n")
