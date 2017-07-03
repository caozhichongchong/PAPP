#set working directory
setwd('yourdir')

#import table of completeness
temp = list.files(pattern="*blast.txt")
result<-list()

#import Cutoff for Core-genome from step 1 Pancutoff.R
Cutoff=Cutoffcore
#----------read blast result and combine into matrix a-----------------------
#read the blast result
for (m in 1:length(temp))
{
  result[[m]]=read.table(temp[m],header=T)
  result[[m]]=as.matrix(result[[m]])
##attention! adjust the character '-' and '_' to simplify the ORF name into the genome name
  for (i in 1:nrow(result[[m]]))
  {
    c<-strsplit(result[[m]][i,2],'-',fixed=TRUE)
    c<-strsplit(c[[1]],'_',fixed=TRUE)
    result[[m]][i,2]<-as.character(c[[1]][1])
  }
}
#combine into matrix a
a<-matrix(0,nrow=(nrow(result[[1]])+nrow(result[[2]])+nrow(result[[3]])),ncol=3)
colnames(a)<-c('Gene&Contig','list_ID','Matrix_ID')
aROW<-0
for (m in 1:(length(temp)))
{
  for (i in 1:nrow(result[[m]]))
  {
    aROW=aROW+1
    a[aROW,1]<-paste(result[[m]][i,1],result[[m]][i,2])
    a[aROW,2]<-m
    a[aROW,3]<-i
  }  
}
##remove duplicate
Unique<-subset(a, !duplicated('Gene&Contig'))
Unique<-Unique[!duplicated(Unique[,1]),]

##store blast unique results
Blastnew=matrix(0,nrow=nrow(Unique),ncol=ncol(result[[1]]))
for (Row in 1:nrow(Unique))
{
  m<-as.numeric(Unique[Row,2])
  i<-as.numeric(Unique[Row,3])
  Blastnew[Row,]=as.matrix(result[[m]][i,])
}

colnames(Blastnew)=c('Query_id','Subject_id','identity','alignment_length'
                     ,'mismatches','gap_openings','q._start','q._end','s._start','s._end','e_value','bit_score','Length_gene','ORF_length','Ratio')
#output blast result
write.table(Blastnew,file='BlastUnique.txt',row.names=F) ###save blast results

#----------------------------Calculate gene frequency------------------
Nameuni<-unique(Blastnew[,1])
Nameuni=unique(Nameuni)
Nameuni=as.matrix(Nameuni)
#gene occurence can be infered from blast results >>> how many different genomes
Count<-matrix(0,nrow=length(Nameuni),ncol=3)
#Nameuni<-levels(Nameuni)[Nameuni]
Count[,1]=Nameuni
unique(Blastnew[,2])
for (i in 1:length(Nameuni))
{
  Count[i,2]=length(which(Blastnew[,1]==Nameuni[i])) ###count occurence
  Check=Blastnew[which(Blastnew[,1]==Nameuni[i]),] ###check result
  if(length(Check)>ncol(Blastnew)){
    Count[i,3]=paste(c(unique(Check[,2])),collapse=" ")
  }else
    Count[i,3]=paste(c(unique(Check[2])),collapse=" ")
}
#output gene frequency
write.table(Count,file='Gene_frequency.txt',row.names = F)
#---------------------Classify pan-genome-----------------------
#select core genes
Core=Count[which(as.numeric(Count[,2])>=Cutoff),]
#dispensable-genome and strain-specific-genome classification
Unique=as.matrix(Count[which(as.numeric(Count[,2])==1),])
Dispen=as.matrix(Count[which(as.numeric(Count[,2]) %in% c(seq(2,Cutoff-1,1))),])

#output pan-genome classified into three subsets
write.table(Core,file='Core_genome.txt', row.names = F, quote = T)
write.table(Dispen,file='Dispensable_genome.txt', row.names = F, quote = T)
write.table(Unique,file='Strain_Specific_genome.txt', row.names = F, quote = T)
