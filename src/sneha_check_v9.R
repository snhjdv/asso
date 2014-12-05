#We run tests and also check covaraice calc by hand derivd formuls: cov(A)=n/(n-1)cov(T)cov(Z)
#We also calculate lambda for each individual and then take mean.
# not used the gam command


#########Start the clock!
ptm <- proc.time()
require(MASS)
require(SKAT)
require(fda)
require(mgcv)
source(file="trait_func.r")

ped<-read.table("ped.ped",header=T,sep="\t");
info<-read.table("info.info",head=T,sep="\t");
#seg.pos<-read

frq<-c(0,.005,.01,.05,.1,.5)

for(fre in 1: length(frq))
{
###function parameters###
nS<-1000;
nN<-500
ntrait<-3
traits<-c("normal","normal","normal");
Beta.miu=0; Beta.sd=10; miu=0;
cmn.frq<-0;   #MAF freq
func.frq<-frq[fre];
nSeq<-3e4;    #length of the sequence


#####initiating variables for the code###

traitschar<-NULL;
for(i1 in 1:ntrait)
{
  if(traits[i1]=="normal"| traits[i1]=="cauchy" | traits[i1]=="T")
    traitschar[i1]<-"C";
  if(traits[i1]=="binary"| traits[i1]=="poisson") traitschar[i1]<-"D";
  
}


cmn = c(cmn.frq, 1)
idx.cmn <- (info$maf>cmn[1] & info$maf<cmn[2])
ustat<-NULL;
cov_A<-NULL;
cov_A_b<-NULL;
pval<-NULL
pval_B<-NULL;
deci_S<-NULL;
pval_S<-NULL;
Results<-data.frame(NULL);

for(i in 1:nS)
{
  ## sliding window######
  seg.pos <- runif(1,0,(max(info$pos)-nSeq));
  idx<-( info$pos>seg.pos & info$pos<seg.pos+nSeq );
  idx<- idx & idx.cmn;
  smp<-sample(1:nrow(ped),nN);
  geno<-ped[smp,idx];
  vrt<-apply(geno,2,variant);#false means everything is same and hence taken out
  indx <-which(idx=="TRUE")[vrt];
  geno<-as.matrix(geno[,vrt]); 
  
  
  ####generating phenotype####
  Tn<-NULL;
  idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq); # index functional alleles
  geno.true <- as.matrix(geno[,idx.true]);
  for(i1 in 1:ntrait)
  {
    Tn<-cbind(Tn,trait.simu(geno.true,Beta.miu,Beta.sd,miu,cov,dist=traits[i1]));
  }
  
  
  #########Flipping genotype data###########
#   geno.flipped = matrix(geno[,1], ncol=1)
#   
#   k = dim(geno)[2]-1
#   for(i in 1:k){
#     geno.flipped = cbind(geno.flipped, flip.check(geno.flipped[,i], geno[,i+1]))
#   }
  
  # ------------------------------------------------------------------------
  # Analysis based on Cubic-Splines 
  # ------------------------------------------------------------------------
  
  shouldContinue=TRUE;
  pos <-info$pos[indx] # physical position of all variants
  pos <- (pos-pos[1])/(pos[length(pos)]-pos[1])
  pos2 = seq(min(pos), max(pos), length=2*length(pos))
  pos3<-seq(0,1,0.001)
  shouldContinue=TRUE;
#   g=numeric(0);
#   for(i3 in 1:nN){
#     G <- gam(geno[i3,]~s(pos,fx=FALSE, k=-1,bs='cr')) ###cr is for cubic regression
#     g = cbind(g, G$fitted.values)
#   }

  
  if(shouldContinue==TRUE)
  {
  lambda<-rep(0,nN)
   for(j in 1:nN)
   {
     genospline<-try(smooth.spline(pos,geno[j,],nknots=length(pos)))
     if("try-error" %in% class(genospline)){
    lambda[j]<-NA}
     else {lambda[j]<-genospline$lambda}
   }
  }
  
  lambda_all<-mean(lambda,na.rm=TRUE)
  Knots = pos
  norder = 4
  nbasis=length(Knots) + norder - 2
  ybasis=create.bspline.basis(range(Knots), nbasis, norder, Knots)
  Lfdobj= 2
  yfdPar = fdPar(ybasis, Lfdobj, lambda=lambda_all)
  yfd = smooth.basis(pos,t(geno),yfdPar)$fd
  yij.f = eval.fd(pos3,yfd)
  z<-NULL;
  z<-rbind(z,apply(yij.f*.001,MARGIN=2,sum))
  
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
  
  Astat<-apply(A,2,sum)
  # Ustat<-(2/sqrt(nN))*Astat;
  
  ################# Covariance of Astat############
  va<-as.numeric(var(t(z)))
  cov_Astat<-(nN/(nN-1))*(cov(Tn)*va);
  cov_A<-rbind(cov_A,cov_Astat)
  
  rankA<-rankMatrix(cov_Astat)[1];
  Astattr<-as.matrix(Astat,nrow=ntrait)
  chisqst<-(Astat%*%ginv( cov_Astat)%*%Astattr);
  pval[i]<-pchisq( chisqst,rankA, ncp = 0, lower.tail = FALSE, log.p = FALSE);
  
  #----------------------------------------------------------------------------
  ######################Burden technique#############################
  #----------------------------------------------------------------------------
  mialfreq<-info$maf[idx];
  mialfreq<- mialfreq[vrt];
  wt1<-1/(sqrt(mialfreq*(1-mialfreq)))
  zb<-NULL;
  Ab<-NULL;
  for( j in 1:nN)
  {
    
    zb[j]<-geno[j,]%*%wt1;
  }
  
  varzb<-var((zb));
  for( i1 in 1:nN)
  {
    Ab<-rbind(Ab,(sqrt(nN)/(nN-1))*( zb[i1]*ubar[i1,]));
  }
  
  Abstat<-apply(Ab,2,sum);
#   Ubstat<-(2/sqrt(nN))*Abstat;
#   cov_Ub<-matrix(rep(0,ntrait*ntrait),nrow=ntrait,ncol=ntrait);
#   for( i1 in 1:nN)
#   {
#     mat<-(ubar[i1,]%*%t(ubar[i1,]))
#     coeffb<-(4/(nN-1)**2)*(var(zb));
#     cov_Ub<-cov_Ub+mat*as.numeric(coeffb) 
#   }
  
  va_b<-as.numeric(var(zb))
  cov_Astat_b<-(nN/(nN-1))*(cov(Tn)*va_b);
  cov_A_b<-rbind(cov_A_b,cov_Astat_b)

  rankUb<-rankMatrix(cov_Astat_b);
  Abstattr<-as.matrix( Abstat,nrow=ntrait)
  chisqst_b<-( Abstat%*%ginv(cov_Astat_b)%*%Abstattr);
  pval_B[i]<-pchisq( chisqst_b,rankUb, ncp = 0, lower.tail = FALSE, log.p = FALSE);
  
  
  
  #-----------------------------------------------------------------
  ###################SKAT#################################
  #-----------------------------------------------------------------
  
  
  pval_S<-NULL;
  for(i1 in 1:ntrait)
  {
    obj <- SKAT_Null_Model( Tn[,i1] ~ NULL, out_type = traitschar[i1] )
    pval_S[i1]<-Results[i,paste("SKAT",i1)] <- as.numeric(try(SKAT(geno,obj)$p.value))
    
  }
  
  
  deci_SKAT<-function(vec)
  {
    for(i1 in 1:ntrait)
    {
      if( as.numeric(vec[i1])< (.05/ntrait)) {dec=1;}
      if( as.numeric(vec[i1])>= (.05/ntrait)) {dec=0;}
      if(dec==1) {break;}
    }
    return(dec); 
  }
  
  deci_S[i]<-deci_SKAT(pval_S)
  print(i);
}

cpvalues1<-data.frame(Results);
cpvalues2<-data.frame(pval,pval_B);
cpvalues<-data.frame(cpvalues1,cpvalues2)

######################Rejection calc  ( reject null when true)############
#Bonflevel<-.05/NTrait;
count_Bon<-function(pvec)
{
  return(length(which(pvec<.05/ntrait))/nS);
}

count<-function(pvec)
{
  return(length(which(pvec<.05))/nS);
}

rejectrate1<-apply(cpvalues1,2,count_Bon);
rejectrate2<-apply(cpvalues2,2,count);
rejectrate3<-sum(deci_S)/nS;
rejectrate<-append(rejectrate1,rejectrate3)
rejectrate<-append(rejectrate,rejectrate2);
resmat<-as.matrix(rejectrate,nrow=1);
resmat2 <- matrix(resmat, ncol = length(rejectrate), dimnames = NULL)
paravalues<-c(func.frq,nS,nN,ntrait,paste(traits,collapse=":"),Beta.miu,Beta.sd,miu,cmn.frq)
summary<-c(paravalues,resmat2)
write(summary,file=paste("smp",nN,"func.frq",func.frq,"BBBreport_v9.csv"),append=TRUE,sep=",",ncolumns=length(summary));
#write.table(cov_A,file=paste("func.frq",func.frq,"calc_covar_v9.csv"),sep=",")
# write(cpvalues,file=paste("pval_","smp",nN,".csv"),append=TRUE)

}

################## Stop the clock
proc.time() - ptm;



