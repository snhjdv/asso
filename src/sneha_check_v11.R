#### this is a new stat where we consider square of distance

ptm <- proc.time()
require(MASS)
require(SKAT)
require(fda)
require(mgcv)
require(doParallel)
require(foreach)
require(mvShapiroTest)
source(file="src/trait_func_ori.R")

ped<-read.table("dat/ped.ped",header=T,sep="\t");
info<-read.table("dat/info.info",head=T,sep="\t");

frq<-c(0)

for(fre in 1: length(frq))
{
	###function parameters###
	nS<-1000L;
	nN<-500L;
	ntrait<-1
	traits<-c("normal");
	Beta.miu=0; Beta.sd=10; miu=0;
	cmn.frq<-0;	 #MAF freq
	func.frq<-frq[fre];
	nSeq<-3e4;		#length of the sequence
	U<-rep(0,ntrait);
	Ustat<-NULL
	#####initiating variables for the code###
	
	traitschar<-NULL;
	for(i1 in 1L:ntrait)
	{
		if(traits[i1]=="normal"| traits[i1]=="cauchy" | traits[i1]=="T")
			traitschar[i1]<-"C";
		if(traits[i1]=="binary"| traits[i1]=="poisson")
			traitschar[i1]<-"D";
	}
	
	cmn = c(cmn.frq, 1)
	idx.cmn <- (info$maf>cmn[1] & info$maf<cmn[2])
	
	
	seg.pos <- runif(1,0,(max(info$pos)-nSeq));
	idx<-( info$pos>seg.pos & info$pos<seg.pos+nSeq );
	idx<- idx & idx.cmn;
	smp<-sample(1:nrow(ped),nN);
	geno<-ped[smp,idx];
	vrt<-apply(geno,2,variant);#false means everything is same and hence taken out
	indx <-which(idx=="TRUE")[vrt];
	geno<-as.matrix(geno[,vrt]); 
	
	
	##########Fitting function#############
	shouldContinue=TRUE;
	pos <-info$pos[indx] # physical position of all variants
	pos <- (pos-pos[1])/(pos[length(pos)]-pos[1])
	pos2 = seq(min(pos), max(pos), length=2*length(pos))
	pos3<-seq(0,1,0.001)
	shouldContinue=TRUE;
	
	idx.knots = seq(2,length(pos)-1, by=2)
	Knots = c(pos[1], pos[idx.knots], pos[length(pos)])
	norder = 4
	nbasis=length(Knots) + norder - 2
	ybasis=create.bspline.basis(range(Knots), nbasis, norder, Knots)
	Lfdobj= 2
	yfdParL = fdPar(ybasis, Lfdobj)
	
	yfdL = smooth.basis(pos,t(geno),yfdParL)$fd
	yij.fL = eval.fd(pos,yfdL)
	z<-NULL;
	z<-rbind(z,apply(yij.f*.001,MARGIN=2,sum))
	registerDoParallel(cl=4L)
	
	out<-foreach(i = 1L:nS, .combine=rbind) %dopar%
	{
		##################generating phenotypes############
		cutoff<-runif(1)
		Tn<-matrix(0,nrow=nN,ncol=ntrait);
		
		idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq); # index functional alleles
		geno.true <- as.matrix(geno[,idx.true]);
		Tn<-NULL;

		for(i1 in 1:ntrait)
		{
			idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq); # index functional alleles
			geno.true <- as.matrix(geno[,idx.true]);
			Tn<-cbind(Tn,trait.simu(geno.true,Beta.miu,Beta.sd,miu,cov,dist=traits[i1]));
		}
		
		for(i1 in 1:nN)
		{
			if(i1<nN)
			{
				for( j1 in (i1+1):nN)
				{
					genediff<-sum((yij.fL[,i1]-yij.fL[,j1])^2)	
					U<-U+(Tn[i1,]-Tn[j1,])^2*genediff				
				}
			}
		}
		U<-2*U/(nN*(nN-1))
		U
	}
	stopImplicitCluster()
}

proc.time() - ptm;


mvShapiro.Test(out)







