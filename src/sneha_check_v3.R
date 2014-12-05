#this code is to check normality with FDA.


#From this code we can see tht our test followrs MNV when func.frq=0.And test of normality is rejected when func.frq=.5.
#We use mshapiro.test
# keep beta.miu=0 bu beta.sd higher like 10.
#tn's are generated as 1's for beta.miu=2,beta.sd=1
require(fda)
source(file="src/trait_func.r")
nN<-300
ntrait<-1
nS<-600;
probvec<-c(.05)
ustat<-NULL;
func.frq<-0.0;
lb.beta=-1; ub.beta=1; miu=0;
control<-rep(0,nN);
traits<-c("binary")
Astat1<-NULL;
ped<-read.table("dat/ped.ped",header=T,sep="\t");
info<-read.table("dat/info.info",head=T,sep="\t");
info<-info[1:400,];
geno<-ped[,1:400];
# ------------------------------------------------------------------------
# Analysis based on Cubis-Splines and SMALL number of basis functions
# ------------------------------------------------------------------------
vrt<-apply(geno,2,variant); #false means everything is same and hence taken out
geno<-as.matrix(geno[,vrt]);
shouldContinue=TRUE;
# indx <-which(idx=="TRUE")[vrt] #index all variants  
pos <-info$pos[vrt] # physical position of all variants
pos <- (pos-pos[1])/(pos[length(pos)]-pos[1])
pos2 = seq(min(pos), max(pos), length=2*length(pos))

bbasis <- max(floor(length(pos)*0.1),4) # number of basis functions
basis <- create.bspline.basis(norder=4,nbasis=bbasis) 
yfdParS <- fdPar(basis, Lfdobj=2)
yfdS <- smooth.basis(pos,t(geno),yfdParS)$fd
pos3<-seq(0,1,by=.0001)
yij.fS = eval.fd(pos3,yfdS)
z<-NULL;
z<-rbind(z,apply(yij.fS*.0001,MARGIN=2,sum))

for(i in 1:nS)
{
	Tn<-matrix(0,nrow=nN,ncol=ntrait);
  idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq); # index functional alleles
  geno.true <- as.matrix(geno[,idx.true]);
  for(i1 in 1:ntrait)
  {
    
  	repeat
  	{
  		#print(i1);
  		idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq); # index functional alleles
  		geno.true <- as.matrix(geno[,idx.true]);
  		Tn[,i1]<-trait.simu(geno.true,lb.beta,ub.beta,miu,NULL,traits[i1],control);  
  		if(length(which(Tn==0))>(nN/10) && length(which(Tn==1))>(nN/10))
  			break;
  	}
  }
  
  u=matrix(data=0,nrow=nN,ncol=nN);
  ubar<-matrix(data=0,nrow=nN,ncol=ntrait)
  A<-NULL;
  
  for(i1 in 1:nN)
  {
    u<-NULL;
    for(j1 in 1:nN) 
    {
      u<-rbind(u,Tn[i1,]-Tn[j1,]);      
    }
    ubar[i1,]<-apply(u,2,mean)
  }
  
  for(i1 in 1:nN)
  {
    A<-rbind(A,(sqrt(nN)/(nN-1))*( z[i1]*ubar[i1,]));
    
  }
  
  Astat1<-rbind(Astat1,apply(A,2,sum))
  #ustat<-(2/sqrt(nN))*Astat1;
  print(i);
}
shapiro.test(t(Astat1))$p.value
