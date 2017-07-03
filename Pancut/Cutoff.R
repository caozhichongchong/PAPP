library("combinat")
Cutoff<-function(Com,Core_FN_cut=0.01,Dis_FN_cut=0.01,Matrix=T){
  #-------------------------Function-------------------------------
  FNcalculation<-function(Com,N,temp,FNlast){
      for (j in 1:ncol(temp))
      {
        Core_FNadd<-1#total probability of independant events
        Core_FNtemp<-matrix(Com[1:N,2],nrow=1,ncol=N)#matrix for each genome probability
        #incompleteness of one genome set in temp in Core_FNtemp matrix
        #probability of genomes not in temp is not changed
        if(nrow(temp)>1)
          for(k in 1:nrow(temp))
            Core_FNtemp[1,temp[k,j]]<-1-Core_FNtemp[1,temp[k,j]]
        else
          Core_FNtemp[1,temp[1,j]]<-1-Core_FNtemp[1,temp[1,j]]
        #total probability of incompleteness
        for(num in 1:N)
          Core_FNadd<-Core_FNadd*Core_FNtemp[1,num]
        #add up for all genome sets in temp
        FNlast=FNlast+Core_FNadd
      }
    return(FNlast)
  }
  #-------------------------Definition of variables-----------------------------
  Result=list()
  Cutoff=0
  N<-nrow(Com)
  #-----------------------------Core FN rate-----------------------------------
  { Core_FN<-matrix(1,nrow=N,ncol=2)
    ##Cutoff=N
    Core_FN_N<-1
    for(j in 1:N)
    {
      Core_FN[j,1]<-N-j+1
      Core_FN[1,2]<-Core_FN[1,2]*Com[j,2]
      Core_FN_N<-Core_FN_N*(1-Com[j,2])
    }
    Core_FN[N,2]<-1-Core_FN_N
    ##Core_FN for Cutoff from 2~N
    if(N>2)
    {
      for (i in 1:(N-2))
      { 
        Core_FN[i+1,2]=Core_FN[i,2]
        temp<-combn(N, i)
        Core_FN[i+1,2]=FNcalculation(Com,N,temp,Core_FN[i+1,2])
      }
    }
    Core_FN[,2]=1-Core_FN[,2]
    ##Cutoff1 for Core_FN_cut
    if(length(Core_FN[which(Core_FN[,2]<=Core_FN_cut),])>2){
        Cutoff1<-Core_FN[which(Core_FN[,2]<=Core_FN_cut),][1,1]
    }else if(length(Core_FN[which(Core_FN[,2]<=Core_FN_cut),])==2){
        Cutoff1<-Core_FN[which(Core_FN[,2]>=Core_FN_cut),][1]
    }else if(length(Core_FN[which(Core_FN[,2]<=Core_FN_cut),])==2)
      Cutoff1<-Core_FN[which(Core_FN[,2]<=Core_FN_cut),][1]
  }
  #-----------------------------Dispensable FN rate-----------------------------------
  { Dis_FN<-matrix(0,nrow=N-2,ncol=2)
    ##Cutoff=N
    for(j in 1:(N-2))
    {
      Dis_FN[j,1]<-N-j
    }
    ##Dis_FN for Cutoff from 1~N
      for (i in 1:(N-2))
      { 
        if(i>1)
          Dis_FN[i,2]=Dis_FN[i-1,2]
        temp<-combn(N, Dis_FN[i,1])
        for(j in 1:ncol(temp))
        {
          temp2=combn(Dis_FN[i,1],(Dis_FN[i,1]-1))
          Comtemp=Com[c(temp[,j]),]
          Dis_FN[i,2]=FNcalculation(Comtemp,Dis_FN[i,1],temp2,Dis_FN[i,2])
        }
        Dis_FN[i,2]=Dis_FN[i,2]/ncol(combn(N, Dis_FN[i,1]))
      }
    ##Cutoff2 for Dis_FN_cut
    if(length(Dis_FN[which(Dis_FN[,2]<=Dis_FN_cut),])>2){
      Cutoff2<-Dis_FN[which(Dis_FN[,2]<=Dis_FN_cut),][1,1]
    }else if(length(Dis_FN[which(Dis_FN[,2]<=Dis_FN_cut),])==2){
      Cutoff2<-Dis_FN[which(Dis_FN[,2]<=Dis_FN_cut),][1]
    }else if(length(Dis_FN[which(Dis_FN[,2]<=Dis_FN_cut),])==2)
      Cutoff2<-Dis_FN[which(Dis_FN[,2]<=Dis_FN_cut),][1]
    ##Cutoff of core for both creteria of FN
    Result[[1]]=min(Cutoff1, Cutoff2)
  }
  #result storing for FN matrix
  if(Matrix==T)
  {
    Result[[2]]=Core_FN
    Result[[3]]=Dis_FN
  }
  return(Result)
}