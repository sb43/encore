
#### This scripts takes raw counts for sgRNA, normalizes them (using the sum per sample), 
#### calculates fold changes and scales them based on positives and negatives.
#### Concerning the ratios, they are calculated using the plasmid as well as parental lines (Cas9 negs, they will be considered separated) 
#### Information for plasmid and parentals should be provided in the form of a table with two columns: sample names in the first one, 
#### name for the control in the second (for example 'plasmid', 'parental1'...). Biological replicas should have the same name (to average them).

#### The following arguments are needed to run the script:
#### 1.-Path and name of the file with the raw counts
#### 2.-Path and name of the file with information regarding controls described above
#### 3.-Name of the library (will be used to create the name for the output files)
#### 4.-Path for the output

#### Regarding the output, several matrixes will be created:
#### 1.- Norm counts
#### 2.- logFC, one per control plus the average of the parental lines (if requested)
#### 3.- scaled logFC (same number as logFC files)
############################################################

# Check if enough arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) < 4) {
  stop("Insufficient arguments. Please provide.")
}

# Get input and output file paths from command line arguments

input_file <- commandArgs(trailingOnly = TRUE)[1]
ctrl_file <- commandArgs(trailingOnly = TRUE)[2]
lib_name <- commandArgs(trailingOnly = TRUE)[3]
output_path <- commandArgs(trailingOnly = TRUE)[4]


# Check if the input file  and the lib file exists
if (!file.exists(input_file)) {
  stop("Input file does not exist.")
}
if (!file.exists(ctrl_file)) {
  stop("Ctrl file does not exist.")
}

####Load the counts and the control information 

counts=as.matrix(read.delim(input_file))
ctrl=as.matrix(read.delim(ctrl_file))

####New matrix containing the ID and other information regarding the guides

IDs=counts[,c("ID","Note","MyNote","sgRNA1_ID","sgRNA2_ID","Gene1","Gene2")]

####Extracting the columns with the controls in to a separate numeric mat

ctrl_mat=counts[,colnames(counts)%in%ctrl[,1]]
ctrl_mat=apply(ctrl_mat,2,as.numeric)
rownames(ctrl_mat)=IDs[,'ID']

####Keeping in 'counts' the numeric columns that are not controls

counts=counts[,9:ncol(counts)]
counts=counts[,!colnames(counts)%in%colnames(ctrl_mat)]
counts=apply(counts,2,as.numeric)
rownames(counts)=IDs[,'ID']

####Normalization of 'counts'

for (i in 1:ncol(counts)){
  counts[,i]=counts[,i]*10000000/sum(counts[,i]) 
}

####Saving  the matrix with norm counts, we will add extra columns in the following loops with the counts for each control

norm_mat_final=counts

################################################################################################################################################
################################################################################################################################################

###### Generating logfoldchanges and scaling for each control 

##Vector containing each type of control (plasmid, parental 1, parental 2...)
controls=unique(ctrl[,2])

for (j in 1:length(controls)){
  
  ##norm of control (if there are biological replicas)
  y=apply(as.matrix(ctrl_mat[,ctrl[ctrl[,2]==controls[j],1]]),1,mean)
  y=y*10000000/sum(y)
  
  #adding column tor norm matrix (for saving after the loop)
  norm_mat_final=cbind(y,norm_mat_final)
  colnames(norm_mat_final)[1]=controls[j]
  
  ##Matrix for logFC & scaling
  mat_lfc=matrix(0,nrow(counts),ncol(counts))
  colnames(mat_lfc)=colnames(counts)
  rownames(mat_lfc)=rownames(counts)
  
  for (i in 1:ncol(mat_lfc)){
    mat_lfc[,i]=log2((counts[,i]+1)/(y+1))
  }
  
  ###Saving matrix (we will reuse it for the scaling)
  write.table(cbind(IDs,mat_lfc), paste(output_path,lib_name,'_log2FC_',controls[j],'_sgRNA.txt',sep=''),quote = F,sep = "\t")
  
  ##Scaling
  
  negcontrol=rep(0,ncol(mat_lfc))
  poscontrol=rep(0,ncol(mat_lfc))
  
  for (i in 1:ncol(mat_lfc)){
    negcontrol[i]=median(mat_lfc[IDs[,'Note']=='NegativeControls',i])
    poscontrol[i]=median(mat_lfc[IDs[,'Note']=='PositiveControls',i])
  }
  
  negpos = negcontrol-poscontrol
  
  for (i in 1:ncol(mat_lfc)){
    mat_lfc[,i]=(mat_lfc[,i]-negcontrol[i])/negpos[i]
  }
  
  ###Saving matrix 
  write.table(cbind(IDs,mat_lfc), paste(output_path,lib_name,'_log2FC_scaled_',controls[j],'_sgRNA.txt',sep=''),quote = F,sep = "\t")
  
}

#####Saving  the norm counts 

write.table(cbind(IDs,norm_mat_final), paste(output_path,lib_name,'_normcount_sgRNA.txt',sep=''), 
            row.names=F,quote = F,sep = "\t")



