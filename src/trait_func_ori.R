###This was modified to add  t-distribution
### this was modified to control variation between phenotypes

variant <- function(x)
{
  length(unique(x))>1
}

trait.simu <- function(geno,mean,sd,miu,
                       cov=NULL, dist=c("normal","cauchy","poisson","binary","T"))
{
  dist <- match.arg(dist);
  
  #beta<-rnorm(ncol(geno),mean,sd)*abs(log10(maf));
  #set.seed(randseed);
  beta<-rnorm(ncol(geno),mean,sd);
  alpha<-c(0,1);
  
  if(is.null(cov)){
    y.tmp<-geno%*%beta
  }else{
    y.tmp<-geno%*%beta;
  }
  if(dist=="normal"){
    y<-miu+y.tmp+rnorm(nrow(geno),0,1)
  }else if(dist=="cauchy"){
    y<-miu+y.tmp;
    y<-rcauchy(length(y),y,)
  }else if(dist=="poisson"){
    y<-abs(-miu-y.tmp)###exp(-miu-y.tmp);
    #y<-#rzipois(length(y), y, pstr0 =0.6) #zero inflated Poisson
    y<-rpois(length(y), y)
  }else if(dist=="binary"){
    y<-1/(1+exp(-miu-y.tmp));
    # y<-as.numeric(y>runif(nrow(geno)))
    y<-round(y);
  }else if(dist=="T"){
    y<-miu+y.tmp;
    y<-rt(length(y),1,y)
  }
  y
}

trait.simu2<-function(geno,mean,sd,miu,
                      cov=NULL, dist=c("normal","cauchy","poisson","binary"),Ntrait=2)
{
  dist <- match.arg(dist);  
  idx.true<-sample(1:ncol(geno),ncol(geno)*func.frq);
  geno<-as.matrix(geno[,idx.true]);
  
  #beta<-rnorm(ncol(geno),mean,sd)*abs(log10(maf));
  #beta<-rnorm(ncol(geno),mean,sd);
  beta<-runif(ncol(geno),mean-sd,mean+sd);
  #alpha<-c(-1,1);
  
  if(is.null(cov)){
    y.tmp<-geno%*%beta
  }else{
    #alpha<-rep(c(-1,1),ncol(cov)/2);
    y.tmp<-geno%*%beta;
  }
  
  Y<-NULL;
  
  for(i in 1:Ntrait){
    if(dist=="normal"){
      y<-miu+y.tmp+rnorm(nrow(geno),0,1)
    }else if(dist=="cauchy"){
      y<-miu+y.tmp;
      y<-rcauchy(length(y),y,)
    }else if(dist=="poisson"){
      y<-exp(-miu-y.tmp);
      # y<-rzipois(length(y), y, pstr0 =0.6)
      y<-rpois(length(y), y)
    }else if(dist=="binary"){
      y<-1/(1+exp(-miu-y.tmp));
      y<-as.numeric(y>runif(nrow(geno)))
    }
    Y<-cbind(Y,y)
  }
  Y
}


Ucalc<-function(a1,n1)  #chosen asymmetric function is identity
{
  M=NULL;
  for( s1 in 1:n1)
    
  {
    
    U=apply(a1,1,'-',a1[s1,]);
    M=cbind(M,t(U));
  }
  M
}







