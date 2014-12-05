### This is to check if the area under squared stat follows normal distributin


#########Start the clock!
ptm <- proc.time()
require(MASS)
require(SKAT)
require(fda)
require(mgcv)
source(file="src/trait_func.r")
require(mvShapiroTest)
ped<-read.table("dat/ped.ped",header=T,sep="\t");
info<-read.table("dat/info.info",head=T,sep="\t");
#seg.pos<-read
frq<-c(0.1)#, 0.01, 0.05, 0.1)
Astat<-NULL


	###function parameters###
	nS<-150;
	nN<-300;
	ntrait<-1
	traits<-c("binary");
	traits_bin<-rep(4,ntrait);
	lb.beta=-1; ub.beta=1; miu=0;

	cmn.frq<-0;   #MAF freq
	func.frq<-0.5;
	nSeq<-3e4;    #length of the sequence
	
	
	#####initiating variables for the code###
	
	traitschar<-NULL;
	for(i1 in 1:ntrait)
	{
		if(traits[i1]=="normal"| traits[i1]=="cauchy" | traits[i1]=="T" | traits[i1]=="laplace")
		{
			traitschar[i1]<-"C";
			traits_bin[i1]<-0;
		}
		if(traits[i1]=="binary"| traits[i1]=="poisson")
		{
			traitschar[i1]<-"D";
			traits_bin[i1]<-1;
		}
		
	}
	
	
	cmn = c(cmn.frq, 1)
	idx.cmn <- (info$maf>cmn[1] & info$maf<cmn[2])


	######specifying var/cov for phenotype cov control
	#    std1<-5; std2<-5; std3<-5
	#    sig<-matrix(c(1*(std1^2),.5*std1*std2,.6*std1*std3,
	#                .5*std1*std2,1*(std2^2),.2*std2*std3,
	#                .6*std1*std3,.2*std2*std3,1*(std2^3)),nrow=ntrait,ncol=ntrait)
	# 
	

	seg.pos <- runif(1,0,(max(info$pos)-nSeq));
	idx<-( info$pos>seg.pos & info$pos<seg.pos+nSeq );
	idx<- idx & idx.cmn;
	smp<-sample(1:nrow(ped),nN);
	geno<-ped[smp,idx];
	vrt<-apply(geno,2,variant);#false means everything is same and hence taken out
	indx <-which(idx=="TRUE")[vrt];
	geno<-as.matrix(geno[,vrt]); 
	maf<-info$maf

	##########Fitting function#############
	shouldContinue=TRUE;
	pos <-info$pos[indx] # physical position of all variants
	pos <- (pos-pos[1])/(pos[length(pos)]-pos[1])
	pos2 = seq(min(pos), max(pos), length=2*length(pos))
	intrvllen<-0.0001
	pos3<-seq(0,1,intrvllen)
	shouldContinue=TRUE;
	
	idx.knots = seq(2,length(pos)-1, by=2)
	Knots = c(pos[1], pos[idx.knots], pos[length(pos)])
	norder = 4
	nbasis=length(Knots) + norder - 2
	ybasis=create.bspline.basis(range(Knots), nbasis, norder, Knots)
	Lfdobj= 2
	yfdParL = fdPar(ybasis, Lfdobj)
	
	yfdL = smooth.basis(pos,t(geno),yfdParL)$fd
	yij.fL = as.matrix(eval.fd(pos3,yfdL))
	yij.flT<-t(yij.fL)

	# pairwise genotype difference 
    genodiff<-matrix(0.0, nN, nN)
	for(i1 in 1L:nN)
	{
		for(j1 in 1L:nN)
		{
			d<-yij.flT[i1,]-yij.flT[j1,]
			genodiff[i1,j1]<-crossprod(d,d);
		}
	}
	genofiff<-genodiff*intrvllen;

for(i in 1:nS)
{
	## sliding window######

	####generating phenotype####
	
	cutoff<-runif(1)
	Tn<-matrix(0,nrow=nN,ncol=ntrait);
	
	idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq); # index functional alleles
	geno.true <- as.matrix(geno[,idx.true]);
	control<-rep(0,nN);
	bin.pos<-which(traits=="binary")
	rem.pos<-which(traits!="binary")
	if(length(bin.pos)>0)
	{
		for( i1 in bin.pos)
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
	}
	
	for(i1 in rem.pos)
	{
		idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq); # index functional alleles
		geno.true <- as.matrix(geno[,idx.true]);
		Tn[,i1]<-trait.simu(geno.true,lb.beta,ub.beta,miu,cov,dist=traits[i1],control); 			
	}
	
	A<-rep(0,ntrait);
	for(i1 in 1L:nN)
	{
		for(j1 in i1:nN)
		{
			A<-A+(Tn[i1,]-Tn[j1,])*genodiff[i1,j1];
		}
	}
	Astat=rbind(Astat,2*A/(nN*(nN-1)))
	# Ustat<-(2/sqrt(nN))*Astat;
	print(i);
}

#mvShapiro.Test(Astat)
shapiro.test(Astat)
hist(Astat)
mean(Astat)
################## Stop the clock
proc.time() - ptm;
