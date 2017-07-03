#set working directory
setwd('yourdir')

#install packages
install.packages("combinat")
library("combinat")

#install function Cutoff
#load Cutoff.R

#import table of completeness
Com=read.table('Example.txt',header=T)
Com=read.table('Your completeness table',header=T)
#set Creteria
Creteria=0.01

#------------------Calculate FN and FP rates for n from 1 to N---------------
N<-nrow(Com)
Result=Cutoff(Com,Creteria,Creteria,T)
#Cutoff for Core-genome under the Creteria
Cutoffcore=Result[[1]]
#Matrix of FN for Core-genome (FP for Dispensable-genome)
Core_FN=Result[[2]]
#Matrix of FN for Dispensable-genome (FP for Strain-specific-genome)
Dis_FN=Result[[3]]

#output FN and FP matrix
write.table(Core_FN, file="Core_FN.txt", sep="\t",col.names =FALSE,row.names=FALSE,quote = FALSE)
write.table(Dis_FN, file="Dis_FN.txt", sep="\t",col.names =FALSE,row.names=FALSE,quote = FALSE)
