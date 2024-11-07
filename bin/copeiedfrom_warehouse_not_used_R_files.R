########

########Vamos a empezar generando una seleccion que incluya el pklasmido y unas cuantas de las muestras que sabemos han llegado hasta el final, 
########Tambien parental lines para probar que esto tira como manda Dios


setwd('/Users/ih6/projects/ENCORE/pipeline/tables_metadata')
###Vamos a rescatar la tabla REAL 

scp  ih6@farm5-login://lustre/scratch126/casm/team215mg/ek11/ENCORE/ENCORE_CURRENT/SAMPLES_ENCORE_ALL_FEB_2023.xlsx .

######Y a crear nuestro propio congelador
#####Primero para todo lo que parece que hay aqui y luego ya para el test

####PARTIMOS DEL EXCEL JUNIO2023



template=as.matrix(read.delim('/Users/ih6/projects/ENCORE/lustre_files/DataFreeze4_Index.txt'))

final_files=as.matrix(read.delim('/Users/ih6/projects/ENCORE/pipeline/tables_metadata/SAMPLES_ENCORE_ALL_JUNE_2023_VBIDS_FINAL.txt'))
all_files=  as.matrix(read.delim('/Users/ih6/projects/ENCORE/pipeline/tables_metadata/SAMPLES_ENCORE_ALL_JUNE_2023_VBIDS_ALL.txt'))
all_files=all_files[!is.na(all_files[,'INDEX']),]

sum(colnames(template)%in%colnames(all_files))

###cut to same columns

all_files=all_files[all_files[,'CELL_LINE_NAME']%in%c('HT-29'),colnames(all_files)%in%colnames(template)]
all_files[,'INDEX']=gsub(' ','',all_files[,'INDEX'])
all_files[,'INDEX']=paste('N',all_files[,'INDEX'],sep='')
all_files=all_files[all_files[,'LOCATION']=='ENCORE_COLO1',]

all_files[grep('.gz',all_files[,"LANE1.READ1"],invert=T),"LANE1.READ1"]=paste(all_files[grep('.gz',all_files[,"LANE1.READ1"],invert=T),"LANE1.READ1"],'.gz',sep='')
all_files[grep('.gz',all_files[,"LANE2.READ1"],invert=T),"LANE2.READ1"]=paste(all_files[grep('.gz',all_files[,"LANE2.READ1"],invert=T),"LANE2.READ1"],'.gz',sep='')
all_files[grep('.gz',all_files[,"LANE1.READ2"],invert=T),"LANE1.READ2"]=paste(all_files[grep('.gz',all_files[,"LANE1.READ2"],invert=T),"LANE1.READ2"],'.gz',sep='')
all_files[grep('.gz',all_files[,"LANE2.READ2"],invert=T),"LANE2.READ2"]=paste(all_files[grep('.gz',all_files[,"LANE2.READ2"],invert=T),"LANE2.READ2"],'.gz',sep='')

write.table(all_files, "/Users/ih6/projects/ENCORE/pipeline/tables_metadata/DataFreeze4_Index_test.txt", 
            row.names=F, quote = F,sep = "\t")

######Preparamos tambien el txt con el nombre de las files que vamoas a necesitar para el test

filtro=c(all_files[,"LANE1.READ1"],all_files[,"LANE1.READ2"],
         all_files[,"LANE2.READ1"],all_files[,"LANE2.READ2"])

write.table(filtro, "/Users/ih6/projects/ENCORE/pipeline/tables_metadata/file_list.txt", 
            row.names=F,  col.names=F,quote = F,sep = "\t")

write.table(filtro[1:4], "/Users/ih6/projects/ENCORE/pipeline/tables_metadata/file_list_test.txt", 
            row.names=F,  col.names=F,quote = F,sep = "\t")


##########
##########Vamos a ver si esto esta saliendo como debe

counts=as.matrix(read.delim('/Users/ih6/projects/ENCORE/pipeline/tables_metadata/All.counts'))
stats=as.matrix(read.delim('/Users/ih6/projects/ENCORE/pipeline/tables_metadata/All.stats'))

counts=counts[2:nrow(counts),]

head(counts[,1])
emre=as.matrix(read.delim('/Users/ih6/projects/ENCORE/pipeline/tables_metadata/emre/All.counts'))

sum(colnames(counts)%in%colnames(emre))

emre=emre[,colnames(counts)]
counts=counts[2:nrow(counts),]

sum(duplicated(emre[,1]))

rownames(emre)=emre[,1]
rownames(counts)=counts[,1]

emre=emre[rownames(counts),]

plot(as.numeric(counts[,2]),
     as.numeric(emre[,2]))

head(counts[,2])
head(emre[,2])


########################################################################
####Vamos a hacer un pequeno test en esta comporacion 

emre=as.matrix(read.delim('/Users/ih6/projects/ENCORE/pipeline/tables_metadata/emre/All.counts'))

temp=emre[,1]
emre=apply(emre[,2:ncol(emre)],2,as.numeric)
rownames(emre)=temp

temp=counts[,1]
counts=apply(counts[,2:ncol(counts)],2,as.numeric)
rownames(counts)=temp

emre=emre[,colnames(counts)]

cor(as.numeric(counts),
    as.numeric(emre))

####Todo esta bien! 
####Vamos. hacer un pequeno test con la anotacion 

head(as.matrix(read.delim('/Users/ih6/projects/ENCORE/lustre_files/release4/LIBS/COLO1_2G.index.sorted.txt')))
head(as.matrix(read.delim('/Users/ih6/projects/ENCORE/lustre_files/release4/LIBS/COLO1_index.sorted.txt')))
head(as.matrix(read.delim('/Users/ih6/projects/ENCORE/lustre_files/release4/LIBS/COLO1_vectorIDS.sorted.txt')))

test1=as.matrix(read.delim('/Users/ih6/projects/ENCORE/lustre_files/release4/LIBS/COLO1_2G.index.sorted.txt'))
test2=as.matrix(read.delim('/Users/ih6/projects/ENCORE/lustre_files/release4/LIBS/COLO1_index.sorted.txt'))

sum(test1[,1]%in%test2[,1])

################################

####Vamos a abrir el resultado del script 

test=as.matrix(read.delim('/Users/ih6/projects/ENCORE/pipeline/tables_metadata/test_merge.txt'))

################################
################################
###Norm y scaling, files para informacion de controles 

## Anotacion de las cell lines Cas9
cas9=as.matrix(read.delim('/Users/ih6/projects/ENCORE/lustre_files/release4/metadata/Parental_Encore_Lines_2023Aug22.txt'))
pass=as.matrix(read.delim('/Users/ih6/projects/database/cellPassport/model_list_20230801.csv',sep=','))
pass=pass[,1:13]

cas9=cbind('',cas9)
colnames(cas9)[1]='SID'
for (i in 1:nrow(cas9)){
  cas9[i,'SID']=pass[pass[,'model_name']%in%cas9[i,'Cell.Line.Name'],'model_id']
}

IDs_cas9_colo1=paste(cas9[cas9[,6]=='ENCORE.Colo.1','SID'],
                     cas9[cas9[,6]=='ENCORE.Colo.1','CPID'],sep='_')
IDs_cas9_colo2=paste(cas9[cas9[,6]=='ENCORE.Colo.2','SID'],
                     cas9[cas9[,6]=='ENCORE.Colo.2','CPID'],sep='_')
IDs_cas9_colo3=paste(cas9[cas9[,6]=='ENCORE.Colo.3','SID'],
                     cas9[cas9[,6]=='ENCORE.Colo.3','CPID'],sep='_')

####con breast
IDs_cas9_brca1=paste(cas9[cas9[,6]=='ENCORE.Breast.1','SID'],
                     cas9[cas9[,6]=='ENCORE.Breast.1','CPID'],sep='_')
IDs_cas9_brca2=paste(cas9[cas9[,6]=='ENCORE.Breast.2','SID'],
                     cas9[cas9[,6]=='ENCORE.Breast.2','CPID'],sep='_')
IDs_cas9_brca3=paste(cas9[cas9[,6]=='ENCORE.Breast.3','SID'],
                     cas9[cas9[,6]=='ENCORE.Breast.3','CPID'],sep='_')

info_c1=cbind(c("lib.COLO.1",IDs_cas9_colo1),
              c('plasmid','parental1','parental1','parental1','parental2','parental2','parental2'))

info_c2=cbind(c("lib.COLO.2",IDs_cas9_colo2),
              c('plasmid','parental1','parental1','parental1','parental2','parental2','parental2'))

info_c3=cbind(c("lib.COLO.3",IDs_cas9_colo3),
              c('plasmid','parental1','parental1','parental1','parental2','parental2','parental2'))

info_b1=cbind(c("lib.BRCA.1",IDs_cas9_brca1),
              c('plasmid','parental1','parental1','parental1','parental2','parental2','parental2'))

info_b2=cbind(c("lib.BRCA.2",IDs_cas9_brca2),
              c('plasmid','parental1','parental1','parental1','parental2','parental2','parental2'))

info_b3=cbind(c("lib.BRCA.3",IDs_cas9_brca3),
              c('plasmid','parental1','parental1','parental1','parental2','parental2','parental2'))

write.table(info_c1, file = '/Users/ih6/projects/ENCORE/pipeline/testing_Rscripts/colo1_analysis/input/ctrl_info_colo1.txt', 
            row.names=F, quote = F,sep = "\t")












#########