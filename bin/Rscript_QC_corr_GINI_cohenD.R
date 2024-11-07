
##### This is a script to calculate 3 QC related metrics needed for the exclussion of samples:
##### 1.- GINI index 
##### 2.- replica correlation (sgRNA  & gene level)
##### 3.- cohen's D for separation of positive and negative control
#####All the outputs are printed together in a csv table

###Needed arguments to run the script

#### The following arguments are needed to run the script:
#### 1.-Path and name of the file with the runmerged raw counts (the columns will be used as rownames for the output file )
#### 2.-Path and name of the file with the scaled logFC values
#### 3.-Name of the library (will be used to create the name for the output files)
#### 4.-Path for the output

#### Regarding the output, a single csv will be generated containing as many rows as :

####Formulas for GINI index and cohen's D 

gini_index <- function (x, corr = FALSE, na.rm = TRUE){
  if (!na.rm && any(is.na(x))) 
    return(NA_real_)
  x <- as.numeric(na.omit(x))
  n <- length(x)
  x <- sort(x)
  G <- sum(x * 1L:n)
  G <- 2 * G/sum(x) - (n + 1L)
  if (corr) 
    G/(n - 1L)
  else G/n
}

dcohen <- function(set1,set2){
    func=(mean(set1,na.rm=T)-mean(set2,na.rm=T))/sqrt((sd(set1,na.rm=T)^2+sd(set2,na.rm=T)^2)/2)
    return(func)
  }
  
########
# Check if enough arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) < 3) {
  stop("Insufficient arguments. Please provide.")
}

# Get input and output file paths from command line arguments

input_file_raw <- commandArgs(trailingOnly = TRUE)[1]
input_file_lfc <- commandArgs(trailingOnly = TRUE)[2]
lib_name <- commandArgs(trailingOnly = TRUE)[3]
output_path <- commandArgs(trailingOnly = TRUE)[4]

##Get tables 

tab_raw=as.matrix(read.delim(input_file_raw,sep=','))
tab_lfc_rna=as.matrix(read.delim(input_file_lfc,sep=','))

##Formating rawcounts

IDs=tab_raw[,c("ID","Note","MyNote","sgRNA1_ID","sgRNA2_ID","Gene1","Gene2")]
tab_raw=tab_raw[,9:ncol(tab_raw)]
tab_raw=apply(tab_raw,2,as.numeric)
rownames(tab_raw)=IDs[,'ID']

# Output table based on raw counts table

output_colnames=c('GINI_index','cohenD_sgRNA','cohenD_gene')

output_tab=matrix(0,ncol(tab_raw),length(names))
rownames(output_tab)=colnames(tab_raw)
colnames(output_tab)=output_colnames

##Calculating GINI (Counts have to be normalized and log scaled.)
##Also, apply filter if needed at guide level using the lfc table 

for (i in 1:nrow(output_tab)){
temp=tab_raw[,rownames(output_tab)[i]]
temp=temp[!is.na(temp)]
output_tab[i,'GINI_index']=gini_index(temp)
}

## Formating LFC tables for the next steps and calculate gene level (tab_lfc_gen)

IDs=tab_lfc_rna[,c("ID","Note","MyNote","sgRNA1_ID","sgRNA2_ID","Gene1","Gene2")]

tab_lfc_rna=tab_lfc_rna[,!colnames(tab_lfc_rna)%in%colnames(IDs)]
tab_lfc_rna=apply(tab_lfc_rna,2,as.numeric)
rownames(tab_lfc_rna)=IDs[,'ID']

temp=data.frame(tab_lfc_rna,G1=IDs[,'Gene1'],G2=IDs[,'Gene2'],A=IDs[,'Note'])
temp <- aggregate(. ~ G1 + G2 + A , data=temp, median, na.action = NULL, na.rm=T)
tab_lfc_gen=temp[,!colnames(temp)%in%c('G1','G2','A')]
rownames(tab_lfc_gen)=paste(temp[,'A'],temp[,'G1'],temp[,'G2'],sep='_')
rm(temp)

## Calculate the correlation matrix

cor_RNA=t(combn(colnames(tab_lfc_rna),2),0,0)
colnames(cor_RNA)=c('rep_A','rep_B','pearson_corr','spearman_corr')

cor_gen=cor_RNA

for (i in 1:nrow(cor_RNA)){
  
  cor_RNA[i,'pearson_corr']=cor(tab_lfc_rna[,cor_RNA[,'rep_A']],
                                tab_lfc_rna[,cor_RNA[,'rep_B']],
                                use="pairwise.complete.obs",method='pearson')
  cor_RNA[i,'pearson_corr']=cor(tab_lfc_rna[,cor_RNA[,'rep_A']],
                                tab_lfc_rna[,cor_RNA[,'rep_B']],
                                use="pairwise.complete.obs",method='spearman')
  
  cor_gen[i,'pearson_corr']=cor(tab_lfc_gen[,cor_gen[,'rep_A']],
                                tab_lfc_gen[,cor_gen[,'rep_B']],
                                use="pairwise.complete.obs",method='pearson')
  cor_gen[i,'pearson_corr']=cor(tab_lfc_gen[,cor_gen[,'rep_A']],
                                tab_lfc_gen[,cor_gen[,'rep_B']],
                                use="pairwise.complete.obs",method='spearman')
}

#####Write the correlation results to the appropriate location according to the input

saveRDS(paste(output_path,'/corr_table_sgRNA_',lib_name,'.csv',sep=''))






