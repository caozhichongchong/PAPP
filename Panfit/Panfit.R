#set working directory
setwd('yourdir')

#install packages
install.packages("combinat")
install.packages("ggplot2")
library("combinat")
library("ggplot2")
#import table of Genefre
Genefre=read.table('Gene_frequency.txt',header=T)

#import Cutoff for Core-genome from step 1 Pancutoff.R
Cutoff=Cutoffcore
#set Creteria
Creteria=0.01

#install function Cutoff
#load Cutoff.R

#---------------add in the colomn of genome name-----------------
##attention! adjust the character '-' and '_' to simplify the ORF name into the genome name
c<-list()
c<-strsplit(levels(Genefre[,1])[Genefre[,1]],'-',fixed=TRUE)
c2<-matrix(0,nrow=nrow(Genefre),ncol=1)
c3<-strsplit(levels(Genefre[,3])[Genefre[,3]],' ',fixed=TRUE)
Genefre=as.matrix(Genefre)
for (i in 1:nrow(Genefre))
{
  c2[i,1]<-unlist(strsplit(c[[i]][1],'_',fixed=TRUE))[1]
  if(!(c2[i,1] %in% c3[[i]])){
    Genefre[i,3]=paste(Genefre[i,3],c2[i,1])
    Genefre[i,2]=as.integer(Genefre[i,2])+1
  }
}
Genefre=cbind(Genefre,c2)
colnames(Genefre)=c('Gene','Frequency','Sharedgenomes','Genome')
Genefre=data.frame(Genefre)
unique(Genefre$Genome)
write.table(Genefre,'Gene_frequency_adjust.txt',row.names=F)

#-----------------calculate the rarefaction curve-----------
Com=read.table(file='Your completeness table',header = TRUE)
Com=Com[order(Com$Genome),] 
N<-max(as.numeric(Genefre[,2]))
corenum=list()#store core-genome size
uninum=list()#store strain-specific-genome size
pansize=list()#store pan-genome size
for (i in 2:N){
  a<-combn(N, i)
  setnum=1
  uninum[[i-1]]=0
  corenum[[i-1]]=0
  pansize[[i-1]]=0
  if(i<N){
    for (j in 1:ncol(a))
    {
      Comnew=Com[match(a[,j],Com[,1]),]
      Cutoffcore=Cutoff(Com,Creteria,Creteria,F)
      if(!Cutoffcore)
        break
      allgenes=list()
      temppan=0
      for(m in 1:length(a[,j])) ###who is the Number i genomes
      {
        tempcore=0
        tempuni=0
        allgenes[[m]]=Genefre[which(as.numeric(Genefre[,4])==a[,j]),]
        allgenes[[m]]=as.matrix(allgenes[[m]])
        for(n in 1:nrow(allgenes[[m]])){
          #as.integer(unlist(strsplit(allgenes[[m]][1,3],' ')))
          temp=as.integer(strsplit(allgenes[[m]][n,3],' ')[[1]])
          if(length(which(match(a[,j],temp)>0))>=Cutoffcore){
            tempcore=tempcore+1
          }else if(length(which(match(a[,j],temp)>0))==1)
            tempuni=tempuni+1
          temppan=temppan+1/length(which(match(a[,j],temp)>0))
          if(temppan==Inf)
            break
        }
        corenum[[i-1]][setnum]=tempcore
        uninum[[i-1]][setnum]=tempuni
        setnum=setnum+1
      }
      pansize[[i-1]][j]=temppan
    }
  }else
  { j<-1
  Comnew=Com
  Cutoffcore=Cutoff(Com,Creteria,Creteria,F)
  if(!Cutoffcore)
    break
  allgenes=list()
  temppan=0
  for(m in 1:length(a)) ###who is the Number i genomes
  {
    tempcore=0
    tempuni=0
    allgenes[[m]]=Genefre[which(as.numeric(Genefre[,4])==a),]
    allgenes[[m]]=as.matrix(allgenes[[m]])
    for(n in 1:nrow(allgenes[[m]])){
      #as.integer(unlist(strsplit(allgenes[[m]][1,3],' ')))
      temp=as.integer(strsplit(allgenes[[m]][n,3],' ')[[1]])
      if(length(which(match(a,temp)>0))>=Cutoffcore){
        tempcore=tempcore+1
      }else if(length(which(match(a,temp)>0))==1)
        tempuni=tempuni+1
      temppan=temppan+1/length(which(match(a,temp)>0))
    }
    corenum[[i-1]][setnum]=tempcore
    uninum[[i-1]][setnum]=tempuni
    setnum=setnum+1
  }
  pansize[[i-1]][j]=temppan}
}
lapply(pansize, cat, "\n", file="Pan_size.txt", append=TRUE)
lapply(corenum, cat, "\n", file="Core_size.txt", append=TRUE)
lapply(uninum, cat, "\n", file="Strain-specific_size.txt", append=TRUE)

#-----------------draw sampling curves----------------------
data=list()
for(i in 1:3)
  data[[i]]=list()
for(i in 2:N){
  data[[1]][[i-1]]=data.frame(
    x1=rep(i,length(corenum[[i-1]][which(corenum[[i-1]]!=0)])),
    y1= corenum[[i-1]][which(corenum[[i-1]]!=0)])
  data[[2]][[i-1]]=data.frame(
    x2=rep(i,length(uninum[[i-1]][which(uninum[[i-1]]!=0)])),
    y2=uninum[[i-1]][which(uninum[[i-1]]!=0)])
  data[[3]][[i-1]]=data.frame(
    x3=rep(i,length(pansize[[i-1]])),
    y3=pansize[[i-1]]
  )
}
###with fitting
dataall=list()
Mean=list()
for(i in 1:3)
{
  dataall[[i]]=data[[i]][[1]]
  Mean[[i]]=matrix(0,nrow=N-1,ncol=2)
  Mean[[i]][1,2]=mean(data[[i]][[1]][,2])
  Mean[[i]][1,1]=mean(data[[i]][[1]][,1])
  for (n in 3:N){
    dataall[[i]]=rbind(dataall[[i]],data[[i]][[n-1]]) ###merge dataframe
    Mean[[i]][n-1,2]=mean(data[[i]][[n-1]][,2])
    Mean[[i]][n-1,1]=mean(data[[i]][[n-1]][,1])
  }
}

###rarefaction curve for core-genome
i<-1
data<-dataall[[i]] 
Kc<-Mean[[i]][nrow(Mean[[i]]),2]
exponential.model <- glm(log(Mean[[i]][,2]-Kc)~ Mean[[i]][,1])

#~~~~~~~~try to export SD
A<-as.numeric(exponential.model$coefficients[1])
B<-as.numeric(exponential.model$coefficients[2])
x=c(2:(nrow(Mean[[i]])+1)) #2 for core, 3 for uni
y=exp(x*B+A)+Kc
Lm=data.frame(
  x=x,
  y=y
)
g0 <- ggplot() +
  geom_point(data=data, aes(x1,y1),color=cbPalette[1],fill=cbPalette[3],shape=21,size=3,alpha = 0.5,stroke=0.1)+ ###shape=16
  geom_line(data=Lm,aes(x,y),color=cbPalette[8], size = 1)+ #family = gaussian(link="log"),formula = y~exp(-x)
  stat_summary(data=data, aes(x1,y1),fun.y = "mean", color=cbPalette[8], size = 4,shape=15, geom ='point') ###shape=22,stroke
g0 +scale_x_continuous(breaks = 3:N)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  
R=cor(Mean[[i]][2:(N-1),2],y) #Pearson  correlation coefficient R2
#output fitting curve parameters
write.table(data.frame(-1/B,exp(A),Kc,R*R), file='Core_Fitting.txt',sep="\t")

###rarefaction curve for strain-specific-genome
i<-2
data<-dataall[[i]] 
Kc<-Mean[[i]][nrow(Mean[[i]]),2]
exponential.model <- glm(log(Mean[[i]][,2]-Kc)~ Mean[[i]][,1])

#~~~~~~~~try to export SD
A<-as.numeric(exponential.model$coefficients[1])
B<-as.numeric(exponential.model$coefficients[2])
x=c(3:(nrow(Mean[[i]])+1)) #2 for core, 3 for uni
y=exp(x*B+A)+Kc
Lm=data.frame(
  x=x,
  y=y
)
g0 <- ggplot() + 
  geom_point(data=data, aes(x2,y2),color=cbPalette[1],fill=cbPalette[3],shape=21,size=3,alpha = 0.5,stroke=0.1)+ ###shape=16
  geom_line(data=Lm,aes(x,y),color=cbPalette[8], size = 1)+ #family = gaussian(link="log"),formula = y~exp(-x)
  stat_summary(data=data, aes(x2,y2),fun.y = "mean", color=cbPalette[8], size = 4,shape=15, geom ='point') ###shape=22,stroke
g0 +scale_x_continuous(breaks = 3:N)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  
R=cor(Mean[[i]][2:(N-1),2],y) #Pearson  correlation coefficient R2
#output fitting curve parameters
write.table(data.frame(-1/B,exp(A),Kc,R*R), file='Strain-specific_Fitting.txt',sep="\t") 

###pan-fitting ###using strain-specific parameters
Gesize=matrix(0,nrow=N,ncol=1)
for(i in 1:N)
  Gesize[i,1]= nrow(Genefre[which(as.numeric(Genefre[,4])==i),])/Com$Completeness[i]
D=mean(Gesize)
Pan=matrix(0,nrow=N,ncol=2)
for(n in 1:N)
{
  Pan[n,1]=n
  Pan[n,2]=D+Kc*(n-1)+exp(A+B*2)*(1-exp((n-1)*B))/(1-exp(B))
 }
i<-3
data<-dataall[[i]] 
Lm=data.frame(
  x=Pan[,1],
  y=Pan[,2]
)
g0 <- ggplot() + 
  geom_point(data=data, aes(x3,y3),color=cbPalette[1],fill=cbPalette[3],shape=21,size=3,alpha = 0.5,stroke=0.1)+ ###shape=16
  geom_line(data=Lm,aes(x,y),color=cbPalette[8], size = 1)+ #family = gaussian(link="log"),formula = y~exp(-x)
  stat_summary(data=data, aes(x3,y3),fun.y = "mean", color=cbPalette[8], size = 4,shape=15, geom ='point') ###shape=22,stroke
g0 +  scale_x_continuous(breaks = 1:N)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())  
R=cor(Mean[[i]][2:(N-1),2],Pan[3:N,2]) #Pearson  correlation coefficient R2
write.table(data.frame(D,Kc,exp(A+B*2),1/(1-exp(B)),R*R), file='Pan_Fitting.txt',sep="\t") 
