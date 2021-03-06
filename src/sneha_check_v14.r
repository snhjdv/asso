### Check the performance on different GAW genes


#########Start the clock!
ptm <- proc.time()
require(MASS)
require(SKAT)
require(fda)
require(mgcv)
source(file="src/trait_func.r")

ped<-read.table("dat/HIF3A_pedGAW",header=T,sep=",");
ped<-ped[,-1]
info<-read.table("dat/HIF3A_infoGAW",head=T,sep=",");
#seg.pos<-read
path<-'bin/MOST -geno %s -ngeno %d -pheno %s -npheno %d -trait %d -bin %s';
most<-sprintf(path, "ddd", 1000L, "ppp", 500L, 3L, '1 0 0');
frq<-c(0, 0.01, 0.05, 0.1)


for(fre in 1: length(frq))
{
	###function parameters###
	nS<-1000;
	nN<-100;
	ntrait<-3
	traits<-c("binary","binary","binary");
	traits_bin<-rep(4,ntrait);
	lb.beta=-1; ub.beta=1; miu=0;
	cmn.frq<-0;   #MAF freq
	func.frq<-frq[fre];
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
	ustat<-NULL;
	cov_A<-NULL;
	cov_A_b<-NULL;
	pval<-NULL
	pval_B<-NULL;
	deci_S<-NULL;
	pval_S<-NULL;
	deci_Q<-NULL;
	deci_T<-NULL;
	deci_Tpr<-NULL;
	Results<-data.frame(NULL);
	######specifying var/cov for phenotype cov control
	#    std1<-5; std2<-5; std3<-5
	#    sig<-matrix(c(1*(std1^2),.5*std1*std2,.6*std1*std3,
	#                .5*std1*std2,1*(std2^2),.2*std2*std3,
	#                .6*std1*std3,.2*std2*std3,1*(std2^3)),nrow=ntrait,ncol=ntrait)
	# 
	
	phenocor<-matrix(0,nrow=ntrait,ncol=ntrait)
	sigma<-matrix(.7,nrow= ntrait,ncol= ntrait)
	diag(sigma)<-1
	trans<-chol(sigma)
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
		maf<-info$maf
		maf<-maf[idx]
		maf<-maf[vrt]
		wts<-1/sqrt(maf)
		
		####generating phenotype####
		
		cutoff<-runif(1)
		Tn<-matrix(0,nrow=nN,ncol=ntrait);
		
		idx.true <- sample(1:ncol(geno),ncol(geno)*func.frq); # index functional alleles
		geno.true <- as.matrix(geno[,idx.true]);
		# 		control1<-rnorm(nN,0,1)
		# 		control2<-rnorm(nN,0,1)
		# 		control3<-rnorm(nN,0,1)
		# 		control4<-rnorm(nN,0,1)
		#control<-cbind(control1,control2,control3,control4)%*%trans
		control<-rep(0,nN);
		bin.pos<-which(traits=="binary")
		rem.pos<-which(traits!="binary")
		if(length(bin.pos)>0)
		{
			for( i1 in bin.pos)
			{
				repeat
				{
				#	print(i1);
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
		
		phenocor<-phenocor+abs(cor(Tn))
		PID<-1:nN;
		FID<-PID;
		genoty1<-cbind(PID,FID,geno);
		pheno1<-cbind(1:nN,Tn)
		#write.table(pheno1, file="/tmp/phenotype.txt", sep=" ",col.names=c("ID","trait1","trait2","trait3"),row.names=FALSE,quote = F);
	#	write.table(genoty1, file="/tmp/genotype.txt", row.names=FALSE,quote = F);
		
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
		yij.fL = eval.fd(pos,yfdL)
		z<-NULL;
		z<-rbind(z,apply(yij.fL*.0001,MARGIN=2,sum))
		
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
		# print("test done")
		
		
		#--------------------------------------------------------------------
		########################MOST ##################################
		#--------------------------------------------------------------------
# 		traits_bin<-paste(traits_bin, collapse = ' ');
# 		cmd<-paste("bin/MOST -geno /tmp/genotype.txt -ngeno ", nN,
# 				   "-pheno /tmp/phenotype.txt -npheno ",nN," -traits ",ntrait," -covar ",0," -bin", traits_bin);
# 		print(cmd);
# 		system(cmd);  ###RUN IT IN THE TERMINAL
# 		ResDan<-NULL;
# 		ResDan<-read.table("Result.txt", header = T);
# 		pval_Q<-NULL;
# 		pval_T<-NULL;
# 		pval_Tpr<-NULL;
# 		pval_Q<-as.numeric(as.vector(ResDan[-1,3]));
# 		pval_T<-as.numeric(as.vector(ResDan[-1,5]));
# 		pval_Tpr<-as.numeric(as.vector(ResDan[-1,7]));
# 		deciD<-function(vec)
# 		{
# 			for(i1 in 1:length(vec))
# 			{
# 				if( as.numeric(vec[i1])< (.05/length(vec))) {dec=1;}
# 				else {dec=0;}
# 				if(dec==1) {break;}
# 			}
# 			return(dec); 
# 		}
# 		
# 		deci_Q[i]<-deciD(pval_Q);
# 		deci_T[i]<-deciD(pval_T);
# 		deci_Tpr[i]<-deciD(pval_Tpr);
# 		
# 		
# 		
		
		
		#-----------------------------------------------------------------
		###################SKAT#################################
		#-----------------------------------------------------------------
		
		
		pval_S<-NULL;
		for(i1 in 1:ntrait)
		{
			obj <- SKAT_Null_Model( Tn[,i1] ~ NULL, out_type = traitschar[i1] )
			pval_S[i1]<-Results[i,paste("SKAT",i1)] <- as.numeric(try(SKAT(geno,obj)$p.value))
			
		}
		#print("skat done")
		
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
	cpvalues2<-data.frame(pval);
	cpvalues<-data.frame(cpvalues1,cpvalues2)
	#deci_Dan<-data.frame(deci_Q,deci_T,deci_Tpr)
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
	
	rejectrate1<-apply(cpvalues1,2,count_Bon);  #for SKAT
	rejectrate2<-apply(cpvalues2,2,count);    # FOR MY TEST
	rejectrate3<-sum(deci_S)/nS;
	#rejectrate4<-apply(deci_Dan,2,mean); #MOST
	rejectrate<-append(rejectrate1,rejectrate3)#,rejectrate4)
#	rejectrate<-append(rejectrate,rejectrate4)
	rejectrate<-append(rejectrate,rejectrate2);
	resmat<-as.matrix(rejectrate,nrow=1);
	resmat2 <- matrix(resmat, ncol = length(rejectrate), dimnames = NULL)
	paravalues<-c(func.frq,nS,nN,ntrait,paste(traits,collapse=":"),lb.beta,ub.beta,miu,cmn.frq,phenocor[1,2]/nS,dim(geno)[2])
	summary<-c(paravalues,resmat2)
	write(summary,file=paste('out/',"smp",nN,"HIF3A_BBB_v14s.csv", sep=''), append=TRUE,sep=",",ncolumns=length(summary));
	#write.table(phenocor,file=paste('out/',"func.frq",func.frq,"phenocor_BBTT_v10.1.csv"),sep=",")
	# write(cpvalues,file=paste("pval_","smp",nN,".csv"),append=TRUE)
	# print(phenocor)
}

################## Stop the clock
proc.time() - ptm;
