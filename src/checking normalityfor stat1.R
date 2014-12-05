#We run tests and also check covaraice calc by hand derivd formuls: cov(A)=n/(n-1)cov(T)cov(Z)
#We also calculate lambda for each individual and then take mean.
# not used the gam command


#########Start the clock!
ptm <- proc.time()
require(MASS)
require(SKAT)
require(fda)
require(mgcv)
source(file="src/trait_norm_nocontrol.R")

ped<-read.table("dat/ped.ped",header=T,sep="\t");
info<-read.table("dat/info.info",head=T,sep="\t");
#seg.pos<-read




	###function parameters###
	nS<-500;
	nN<-250
	ntrait<-1
	traits<-c("binary");
	Beta.miu=0; Beta.sd=10; miu=0;
	cmn.frq<-0;   #MAF freq
	func.frq<-0.5;
	nSeq<-3e4;    #length of the sequence
	Astat1<-NULL
	
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
	pos3<-seq(0,1,0.0001)
	shouldContinue=TRUE;
	
	idx.knots = seq(2,length(pos)-1, by=2)
	Knots = c(pos[1], pos[idx.knots], pos[length(pos)])
	norder = 4
	nbasis=length(Knots) + norder - 2
	ybasis=create.bspline.basis(range(Knots), nbasis, norder, Knots)
	Lfdobj= 2
	yfdParL = fdPar(ybasis, Lfdobj)
	
	yfdL = smooth.basis(pos,t(geno),yfdParL)$fd
	yij.fL = eval.fd(pos3,yfdL)
	z<-NULL;
	z<-rbind(z,apply(yij.fL*.0001,MARGIN=2,sum))


	
	for(i in 1:nS)
	{
		## sliding window######
	
		
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
		Astat1<-rbind(Astat1,Astat)
		print(i);
	
		
	}
shapiro.test(Astat1)
hist(Astat1)

################## Stop the clock
proc.time() - ptm;

length(Astat1)
