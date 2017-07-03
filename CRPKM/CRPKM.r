#set working directory
setwd('yourdir')

#install packages
install.packages("plyr")
library(plyr)
#
#import RNA and DNA datasets with Name Average.coverage  Reference.length (from assembly)
N=2 #total number of your RNA datasets
RNA1=read.delim('Your RNA1')
RNA2=read.delim('Your RNA2')
DNA1=read.delim('Your DNA')
Reads=read.delim('Reads_number.txt',header=T)
Length=read.delim('DNA_Length.txt',header=T)
Blast=read.delim('blast.txt',header=T,sep='\t')
#Creteria for homology
Identity=50
Hitlength=0.5

#-----------------------------Merge RNA rawreads and calculate RPKM------------------------
#merge all the RNA rawreads into the same dataframe
RNA=merge(data.frame(RNA1),data.frame(RNA2),by='Name')
###calculate RPKM
for(i in 1:nrow(RNA))
{
  RNA$RNA1_rawread[i]=as.numeric(RNA$Average.coverage.x[i])*as.numeric(RNA$Reference.length.x[i])/100
  RNA$RNA2_rawread[i]=as.numeric(RNA$Average.coverage.y[i])*as.numeric(RNA$Reference.length.y[i])/100
}
#output Rawreads
write.table(RNA,'RNA_Rawreads.txt',row.names = F, quote = F)
#calculate RPKM value
for(i in 1:nrow(RNA))
{
  RNA$RNA1_RPKM[i]=as.numeric(RNA$RNA1_rawread[i])/as.numeric(RNA$Reference.length.x[i])/sum(as.numeric(RNA$RNA1_rawread))*1000*1000000
  RNA$RNA2_RPKM[i]=as.numeric(RNA$RNA2_rawread[i])/as.numeric(RNA$Reference.length.y[i])/sum(as.numeric(RNA$RNA2_rawread))*1000*1000000
}
write.table(RNA,'RNA_RPKM.txt',row.names = F, quote = F)
#-------------------------DNA Rawreads----------------------------------------------
for(i in 1:nrow(DNA1))
{
  DNA1$DNA1_RPKM[i]=as.numeric(DNA1$Average.coverage[i])*as.numeric(DNA1$Reference.length[i])/100
}
write.table(DNA1,'DNA_RPKM.txt',row.names = F, quote = F)
#--------------------------MRPKM and CRPKM normalization---------------------------------
DNA=DNA[,-c(2,(ncol(DNA)-1))]
RNA=RNA[,-c(2,(ncol(RNA)-N))]
DNA=merge(DNA,Length[,c(1,3)],by.x='Name')
RNA=merge(RNA,Length,by='Name')
RPKM<-function(Rawreads,Length,Readsnumber){
  RPKM_results=Rawreads
  for(i in 1:nrow(Rawreads))
  {
    for(j in 2:ncol(Rawreads))
    {
      RPKM_results[i,j]=as.numeric(Rawreads[i,j])*10^9/as.numeric(Length[i,2])/as.numeric(Readsnumber[j-1])
    }
  }
  return(RPKM_results)
}
temp1=RNA[,-ncol(RNA)]
temp2=data.frame(Name=RNA$Name,
                 leng1=RNA$Length)
temp3=data.frame(reads1=Reads$RNA1_RPKM,reads2=Reads$RNA2_RPKM)
RNA_RPKM=RPKM(temp1,temp2,temp3)
temp1=DNA[,-ncol(DNA)]
temp2=data.frame(Name=DNA$Name,
                 leng1=DNA$Length)
temp3=data.frame(reads1=Reads$DNA1_RPKM)
DNA_RPKM=RPKM(temp1,temp2,temp3)
RNA_DNA=merge(RNA_RPKM,DNA_RPKM,by='Name',all=T)
MRPKM=RNA_DNA
for(i in 2:(ncol(RNA_DNA)-1))
  MRPKM[,i]=RNA_DNA[,i]/RNA_DNA[,ncol(RNA_DNA)]
#output RNA and DNA RPKM merge table and MRPKM table
write.table(RNA_DNA,'RNA_DNA.txt',row.names = F,sep ='\t',quote=F)
write.table(MRPKM,'MRPKM.txt',row.names = F,sep ='\t',quote=F)
#--------------------------Gene copy number---------------------
#homology
Blast=Blast[which(Blast$identity>=Identity),]
Blast=merge(Blast, Length, by.x='Subject_id',by.y='Name')
Blast=Blast[which(as.numeric(Blast$alignment_length)/as.numeric(Blast$Length)>=Hitlength),]
for(i in 1:nrow(Blast))
{
  ##attention! adjust the character '-' and '_' to simplify the ORF name into the genome name
  Blast$G1[i]=strsplit(toString(Blast$Subject_id[i]),'-')[[1]][1]
  Blast$G2[i]=strsplit(toString(Blast$Query_id[i]),'-')[[1]][1]
}
Blast=Blast[which(Blast$G1==Blast$G2),]
Blast1=count(Blast, c("Subject_id"))
#output Gene copy number per genome
write.table(Blast,'Gene_copy_number_Blast.txt',row.names = F,sep='\t',quote=F)
#--------------------------CRPKM normalization-----------------------
MRPKM=merge(MRPKM,Blast1,by.x='Name',by.y='Query_id',all.x=T)
MRPKM=MRPKM[which(!(is.na(MRPKM$freq))),]
CRPKM=MRPKM[,1:3]
for (i in 2:3)
  CRPKM[,i]=MRPKM[,i]*MRPKM$freq
CRPKM=subset(CRPKM,!duplicated(CRPKM$Name))
MRPKM=subset(MRPKM,!duplicated(MRPKM$Name))
#output CRPKM
write.table(CRPKM,'CRPKM.txt',row.names = F,sep ='\t',quote=F)
write.table(MRPKM,'MRPKM_Copy_Number.txt',row.names = F,sep ='\t',quote=F)
