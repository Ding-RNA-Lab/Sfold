	#passing the parameter in a vector= c("bindsiteFile=value","conserFile=value"(optional), "enrichScoreFile=value","logitCoefFile=value","cdsBegin=value","cdsEnd=value","component=value","type=value") which can have any order; it should not have spaces on both sides of "=";
	
#The function is used for computation of conservation score (when conserFile exists) and relative percentage of siteBegin in each component, then based on the features computed, functional probabilities of each binding site are computed using a nonlinear logistic prediction model.  
#for seed binding site, the conservation score for seed part is computed. 
#for seedless binding site, only the conservation score of whole site is computed.


####### input the parameter #######
	cmd_args = commandArgs(trailingOnly=TRUE); 
	for (x in strsplit(cmd_args,"=")) assign(x[1],x[2]); rm(cmd_args,x);  #convert each paramter in the parameter list to same variable and its value.	

############################ main code ################
	 

	bindsite = read.csv(bindsiteFile,sep="\t",header=T); 
	ind = match(c("SiteID","Transcript","miRNA","UTR_length","Site_range","Seed_range","Seed_kind","region12_17","Ghybrid","Gnucl","Gtotal","SiteAcc","SeedAcc",paste(rep(c("UpsAcc","DwsAcc","UpAU","DwAU"),6),"_",rep(1:6*5,each=4),"nt",sep="")),colnames(bindsite));
	bindsite = bindsite[,ind]; rm(ind); gc();
	
	if (grepl("seedless",type,ignore.case=TRUE)==TRUE) bindsite = bindsite[,!grepl("Seed",colnames(bindsite),ignore.case=T)];
	if (length(ind<-which("Site_range"==colnames(bindsite)))>0) colnames(bindsite)[ind] = "site_position";
	if (length(ind<-which("Seed_range"==colnames(bindsite)))>0) colnames(bindsite)[ind] = "seed_position";

	ind = match("site_position",colnames(bindsite));
	if (is.na(ind)) {print("Warnings: cannot find site position");q(save="no");}
	bsPosit = as.numeric(unlist(strsplit(as.vector(bindsite[,ind]),"-"))); 
	bsPosit = matrix(bsPosit,2);	 
# compute percentage of site start position for current component, if cdsBegin, cdsEnd >0
	cdsBegin = as.integer(cdsBegin); cdsEnd=as.integer(cdsEnd);
if (cdsBegin>0 & cdsEnd>0){
	if (grepl("3pUTR",component,ignore.case=TRUE)) SiteBegin = (bsPosit[1,]-cdsEnd)/(as.vector(bindsite[1,which("UTR_length"==colnames(bindsite))])-cdsEnd);
	if (grepl("CDS",component,ignore.case=TRUE)) SiteBegin = (bsPosit[1,]-cdsBegin+1)/(cdsEnd-cdsBegin+1);
	if (grepl("5pUTR",component,ignore.case=TRUE)) SiteBegin = bsPosit[1,]/(cdsBegin-1);
	bindsite =data.frame(bindsite,SiteBegin); rm(SiteBegin); gc();
	} else if (!exists("conserFile")) {
						SiteBegin=bsPosit[1,]/as.vector(bindsite[1,which("UTR_length"==colnames(bindsite))]);
						bindsite =data.frame(bindsite,SiteBegin); rm(SiteBegin); gc();
						}
#compute conservation score	
	if (exists("conserFile")){
		#search in conservation file for the given gene accession number
		conserScore = readLines(zz<-pipe(paste("cat ",conserFile," | awk '$1~\"",as.vector(bindsite[1,which("Transcript"==colnames(bindsite))]),"\" {print $4,$5}'",sep=""))); close(zz);
		if (length(conserScore)==0) {print(paste(as.vector(bindsite[1,which("Transcript"==colnames(bindsite))]),"cannot be found in trackTable!")); q(save="no");}
		
		conserScore = strsplit(conserScore,"\"");
		cdsInfo = sapply(conserScore,function(x) x[2]);
		conserScore = sapply(conserScore,function(x) x[4]);
		conserScore = strsplit(conserScore,",");
		if (length(conserScore)>1) {
				leng = sapply(conserScore,length);
				i0 = which.min(abs(leng-as.vector(bindsite[1,which("UTR_length"==colnames(bindsite))])));
				conserScore = conserScore[i0];
				cdsInfo = cdsInfo[i0];
				}
		#if the SiteBegin percentage doesnot computed above, compute it using the CDS information from conservation file now.
		if (!is.element("SiteBegin",colnames(bindsite))){
			cdsInfo = as.numeric(unlist(strsplit(cdsInfo,"-")));
			cdsBegin = cdsInfo[1];  cdsEnd = cdsInfo[2]; 
			if (grepl("3pUTR",component,ignore.case=TRUE)) SiteBegin = (bsPosit[1,]-cdsEnd)/(as.vector(bindsite[1,which("UTR_length"==colnames(bindsite))])-cdsEnd);
			if (grepl("CDS",component,ignore.case=TRUE)) SiteBegin = (bsPosit[1,]-cdsBegin+1)/(cdsEnd-cdsBegin+1);
			if (grepl("5pUTR",component,ignore.case=TRUE)) SiteBegin = bsPosit[1,]/(cdsBegin-1);
			bindsite =data.frame(bindsite,SiteBegin); rm(SiteBegin,cdsInfo); gc();
			}
		
		
		conserScore = as.numeric(unlist(conserScore));gc();
		Num= nrow(bsPosit);
		site_conser=rep(NA,Num);  

		#compute seed position
		if (grepl("seedless",type,ignore.case=TRUE)==FALSE){
			ind = match("seed_position",colnames(bindsite));
			if (is.na(ind)) {print("Warnings: cannot find seed position");q(save="no");}
			seed_position = as.vector(bindsite[,ind]);
			seed_position = matrix(as.numeric(unlist(strsplit(seed_position,"-"))),2);
			seed_conser=rep(NA,Num); offseed_conser=rep(NA,Num)
			}; 	
		gc();

		#start to compute the conservation score.
		L = length(conserScore);
		ind = which(bsPosit[2,]<=L);   
		if (length(ind)>0) {
			site_conser[ind] = apply(as.matrix(bsPosit[,ind]),2,function(x,score) sum(score[x[1]:x[2]]),conserScore);
			siteLeng = bsPosit[2,]-bsPosit[1,]+1;
			if (grepl("seedless",type,ignore.case=TRUE)==FALSE){
				seed_conser[ind] = apply(as.matrix(seed_position[,ind]),2,function(x,score) sum(score[x[1]:x[2]]),conserScore);
				seedLeng = seed_position[2,]-seed_position[1,]+1;
				offseed_conser = (site_conser-seed_conser)/(siteLeng-seedLeng);
				seed_conser = seed_conser/seedLeng; rm(seedLeng);
				}
			site_conser = site_conser/siteLeng; rm(siteLeng); gc();
			} 
	
		rm(L,conserScore);  gc(); 
		ind = which(is.na(site_conser)==FALSE);
		bindsite = bindsite[ind,];
		site_conser = site_conser[ind]; 
		#store the conservation information into bindsite
		if (grepl("seedless",type,ignore.case=TRUE)==FALSE){
			rm(seed_position);
			seed_conser = seed_conser[ind]; offseed_conser = offseed_conser[ind];
			bindsite = data.frame(bindsite,site_conser,seed_conser,offseed_conser); 
			rm(site_conser,seed_conser,offseed_conser,ind)
			} else {
					bindsite = data.frame(bindsite,site_conser); rm(site_conser,ind);
				}
	}
	rm(bsPosit); gc()
## end of conservation computing
	
	
# compute logit probability	using nonlinear logistic model
	#find the enriched feature used for model training
	logRegCoef = read.csv(logitCoefFile,header=T,sep="\t",row.names=NULL);
	coefItem = as.vector(logRegCoef[,1]);
	indCoef = which(coefItem=="Intercept")[2]-1;
	selectFeature = unlist(strsplit(coefItem[2:indCoef],"-"));
	selectFeature = selectFeature[seq(1,length(selectFeature),2)];
#since worm project used new script framework, the feature names were also updated. So they need to be converted to old names.
	if (grepl("worm",logitCoefFile,ignore.case=T)) {
		OldFeaName = c("Ghybrid","Gtotal","Gnucl","Seed_kind","region12_17","SiteBegin","SiteAcc","SeedAcc", paste(rep(c("UpsAcc","DwsAcc","UpAU","DwAU"),each=6),rep(paste(1:6*5,"nt",sep=""),4),sep="_"),"site_conser","seed_conser","offseed_conser");
		NewFeaName = c("dG_hybrid","dG_total","dG_nucl","Seed_Type","X3._bp","Site_Location","Site_Access","Seed_Access", paste(rep(c("Upstream_Access","Dwstream_Access","Upstream_AU","Dwstream_AU"),each=6),".",rep(paste(1:6*5,"nt",sep=""),4),".",sep=""),"Site_Consv","Seed_Consv","Offseed_Consv");
		selectFeature = OldFeaName[match(selectFeature,NewFeaName)];
		}
	selectCol = match(selectFeature,colnames(bindsite));
	dataBS = bindsite[,selectCol];
	 
	#compute the enrichment score for seed and 3' complementary part
	enrichScore = readLines(enrichScoreFile);
	if (!is.na(i0<-match("Seed_kind",selectFeature))){
		ind = which(grepl("6mer\t7mer-A1\t7mer-m8\t8mer",enrichScore,ignore.case=T));
		bin = strsplit(enrichScore[ind],"\t",fixed=T)[[1]][2:6]; 
		ratio = strsplit(enrichScore[ind+1],"\t",fixed=T)[[1]][2:6];
		dataBS[,i0] = as.numeric(as.vector(factor(as.vector(dataBS[,i0]),levels=bin,labels=ratio)));
		} 
	if (!is.na(i0<-match("region12_17",selectFeature))){
		ind = which(grepl("Bin      :\t0\t1",enrichScore,ignore.case=T));
		bin = strsplit(enrichScore[ind],"\t",fixed=T)[[1]][2:3];
		ratio = strsplit(enrichScore[ind+1],"\t",fixed=T)[[1]][2:3];
		dataBS[,i0] = as.numeric(as.vector(factor(as.vector(dataBS[,i0]),levels=bin,labels=ratio)));
		} 
##end of enrichment
	
	dataBS = as.matrix(dataBS); gc(); 
	M = nrow(dataBS);
	LogitProb = cbind(1,dataBS)%*%as.matrix(logRegCoef[1:indCoef+indCoef,2:11]);
	#for interaction items
	index = as.vector(logRegCoef[(2*indCoef+1):nrow(logRegCoef),1]); index = matrix(as.integer(unlist(strsplit(gsub("V","",index),"*",fixed=T))),2);
	logRegCoef = as.matrix(logRegCoef[(2*indCoef+1):nrow(logRegCoef),2:11]); N = ncol(index); gc();
	for (i in 1:N) {
		LogitProb = LogitProb + dataBS[,index[1,i]]*dataBS[,index[2,i]]*matrix(rep(logRegCoef[i,],each=M),M);
		gc();
		}
	rm(dataBS,logRegCoef,indCoef); gc();
	LogitProb = exp(LogitProb);
	LogitProb = rowMeans(LogitProb/(1+LogitProb)); 
  #add four columns from right side for "Upstream Access","Dwstream Access","Upstream AU","Dwstream AU" with best window sizes, aimed to show them in StarMir website. If some of them are not enriched, the value is NA. For downloaded, the four columns will be removed.
	UpDw = matrix(NA,M,4); colnames(UpDw) = c("Upstream Access","Dwstream Access","Upstream AU","Dwstream AU");
	selectFeature = selectFeature[grepl("Up",selectFeature) | grepl("Dw",selectFeature)];
	if (length(selectFeature)>0) {
		tmp = matrix(unlist(strsplit(selectFeature,"_",fixed=T)),2);
		ind = match(tmp[1,],c("UpsAcc","DwsAcc","UpAU","DwAU"));
		colnames(UpDw)[ind] = paste(colnames(UpDw)[ind],"(",tmp[2,],")*",sep=""); 
		UpDw[,ind] = as.matrix(bindsite[,match(selectFeature,colnames(bindsite))]); rm(tmp,selectFeature);
		}
	bindsite =data.frame(bindsite,LogitProb,UpDw); 

	
	#rename all features for output
	newName = c("SiteID","Target","miRNA","Target_Len","Site_Position","Seed_Position","Seed_Type","3'_bp","dG_hybrid","dG_nucl","dG_total","Site_Access","Seed_Access",paste(rep(c("Upstream_Access","Dwstream_Access","Upstream_AU","Dwstream_AU"),6),"(",rep(1:6*5,each=4),"nt)",sep=""),"Site_Location","Site_Consv","Seed_Consv","Offseed_Consv","LogitProb",colnames(UpDw));	rm(LogitProb,UpDw); gc();
	if (grepl("seedless",type,ignore.case=TRUE)) newName = newName[!grepl("Seed",newName,ignore.case=T)];
	if (!exists("conserFile")) newName = newName[!grepl("Consv",newName,ignore.case=T)];
	newName[selectCol] = paste(newName[selectCol],"*",sep="");
	colnames(bindsite) = newName;

# save file	
	ind = as.vector(gregexpr("/",bindsiteFile,fixed=T)[[1]]);
	if (ind[1]>0) {ind=ind[length(ind)];bindsiteFile=paste(substr(bindsiteFile,1,ind),"Output-",substring(bindsiteFile,ind+1),sep="")} else bindsiteFile=paste("Output-",bindsiteFile,sep="");
	write.table(bindsite,file=bindsiteFile,sep="\t",col.names=T,row.names=F);
	
	q(save="no");
	

