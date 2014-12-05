#this code is to check kendals tau(with just one genotype)

source(file="src/trait_func.r")
nvariants<-1
nN<-200
ntrait<-3
probvec<-c(.95,.025,.025)
ustat<-NULL;
Beta.miu=2; Beta.sd=1; miu=0;
traits<-c("binary")
Astat1<-NULL;

for(i in 1:500)
{
  genotr<-NULL;
  for(i1 in 1: nvariants)
  {
    
    genotr<-rbind(genotr,sample(c(0,1,2),nN,replace=TRUE,probvec))
  }
geno<-t(as.matrix(genotr,nrow=nN))
z<-geno;
geno.true<-matrix(data=0,nrow=nN);  #To check null hypothesis
Tn<-NULL;
for(i1 in 1:ntrait)
{
  
  Tn<-cbind(Tn,trait.simu(geno.true,Beta.miu,Beta.sd,miu,cov,dist=traits[i1]));
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
