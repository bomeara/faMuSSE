#' @import diversitree
#' @import sfsmisc
#' @import ape
#' @import partitions
#' @import corpcor
#source("V7_UtilityFns.R")



#the following 56 lines are truly stupid. However, otherwise I get an error thrown when i call make.cache.musse
# source('diversitree/R/asr-bisse.R')
# source('diversitree/R/asr-mkn.R')
# source('diversitree/R/asr-musse.R')
# source('diversitree/R/asr.R')
# source('diversitree/R/check.R')
# source('diversitree/R/clade-tree.R')
# source('diversitree/R/constrain.R')
# source('diversitree/R/diversitree-branches.R')
# source('diversitree/R/drop.tip.fixed.R')
# source('diversitree/R/history.R')
# source('diversitree/R/mcmc-norm.R')
# source('diversitree/R/mcmc-slice.R')
# source('diversitree/R/mcmc.R')
# source('diversitree/R/mle-grid.R')
# source('diversitree/R/mle-integer.R')
# source('diversitree/R/mle-mixed.R')
# source('diversitree/R/mle-optimize.R')
# source('diversitree/R/mle-subplexR.R')
# source('diversitree/R/mle-tgp.R')
# source('diversitree/R/mle.R')
# source('diversitree/R/model-bd-ode.R')
# source('diversitree/R/model-bd-split.R')
# source('diversitree/R/model-bd-t.R')
# source('diversitree/R/model-bd.R')
# source('diversitree/R/model-bisse-split.R')
# source('diversitree/R/model-bisse-t.R')
# source('diversitree/R/model-bisse-td.R')
# source('diversitree/R/model-bisse-unresolved.R')
# source('diversitree/R/model-bisse.R')
# source('diversitree/R/model-bm.R')
# source('diversitree/R/model-geosse.R')
# source('diversitree/R/model-mkn.R')
# source('diversitree/R/model-musse-split.R')
# source('diversitree/R/model-musse-t.R')
# source('diversitree/R/model-musse-td.R')
# source('diversitree/R/model-musse.R')
# source('diversitree/R/model-quasse-common.R')
# source('diversitree/R/model-quasse-fftC.R')
# source('diversitree/R/model-quasse-fftR.R')
# source('diversitree/R/model-quasse-mol.R')
# source('diversitree/R/model-quasse-split.R')
# source('diversitree/R/model-quasse.R')
# source('diversitree/R/plot-alt-extra.R')
# source('diversitree/R/plot-alt-util.R')
# source('diversitree/R/plot-alt.R')
# source('diversitree/R/profiles-plot.R')
# source('diversitree/R/simulate-bd.R')
# source('diversitree/R/simulate-bisse.R')
# source('diversitree/R/simulate-musse.R')
# source('diversitree/R/simulate-quasse.R')
# source('diversitree/R/simulation.R')
# source('diversitree/R/split-recycle.R')
# source('diversitree/R/split.R')
# source('diversitree/R/t.R')
# source('diversitree/R/td.R')
# source('diversitree/R/util.R')
# source('diversitree/R/')


#print("Using command script v. Aug 11, 2015")

replaceextralist<-function(i) {
	#print("in replaceextralist")
	#print(paste("i = ",i))
	#print(paste("class i = ",class(i)," dim(i) =",dim(i), " length(i) = ",length(i)))
	assign("extralist",i,envir = .GlobalEnv)
}

make.musse.modifiedWithRootFixedAt1 <- function(tree, states, k, sampling.f=NULL, strict=FALSE,
                       safe=FALSE) { #we turn off strict so that we can still run even if we have just chars 1, 2, 3, 5, 6, 7, 8, for example. And remember that states start with 1, not 0, for musse.
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  branches <- make.branches.musse(k, safe)
  root.p.vector=rep(0,k)
  root.p.vector[1]=1
  ll.musse <- function(pars, condition.surv=TRUE, root=ROOT.GIVEN,
                       root.p=root.p.vector, intermediates=FALSE) {
    if ( length(pars) != k*(k+1) )
      stop(sprintf("Invalid length parameters (expected %d)",
                   k*(k+1)))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    ll.xxsse(pars, cache, initial.conditions.musse, branches,
             condition.surv, root, root.p, intermediates)
  }

  ll <- function(pars, ...) ll.musse(pars, ...)
  class(ll) <- c("musse", "function")
  attr(ll, "k") <- k
  ll
}

#F=focal states
doUnifiedRun<-function(F=F, T=T,D=D,S=6) {
	startTime<-proc.time()
	#first, make the data and tree files
	filename=prepData(P="1_2_3_4_5_6",F=F,T=T,D=D,S=S)
	focalVector=F
	if (length(focalVector==1)) {
		focalVector<-stringToVector(focalVector)
	}
	#next, load the actual data
	data<-read.csv(file=paste(filename,".csv",sep=""),header=TRUE)
	states<-data[,2]
	names(states)<-data[,1]
	phy<-read.tree(file=paste(filename,"tree",".t",sep=""))

	#now start musse setup
	nAngiosperms=250000

	sampling.f<-rep(length(states)/nAngiosperms,2^S)
	results.vector.all<-c()
	lik <- make.musse.modifiedWithRootFixedAt1(tree=phy, states=states, k=2^S, sampling.f=sampling.f)
	#print("trying to assign extralist")
	assign("extralist",list(),envir = .GlobalEnv) #naughty
	#print(paste("first extralist is ",extralist))
	lik.trans <- modify_transitions(lik, type=T, F=F, S=S,extralist=extralist)
	#print(paste("second extralist is ",extralist))
	argnames(lik.trans)
	lik.final <- modify_diversification(lik.trans, type=D, F=F, S=S,extralist=extralist)
	#print(paste("third extralist is ",extralist))
	argnames(lik.final)
	p <- starting.point.musse(phy, 2^S)
	if (length(extralist)>0) {
		p <- starting.point.musse.extra(phy,2^S,argnames=argnames(lik.final)) #this is done as otherwise won't get right starting vector
	}
	fit.final <- find.mle(lik.final,p,method="subplex",hessian=TRUE) #the hessian lets us get standard errors
	fit.final.se<-fit.final
	fit.se<-rep(NA,length(fit.final$par))
	try(fit.se<-sqrt(diag(pseudoinverse(-1*fit.final$hessian)))) #just for extra protection
	#print(fit.se)
	fit.final.se$par<-fit.se
	names(fit.final.se$par)<-names(fit.final$par)
	#print(paste(names(fit.final$par),".se",sep=""))
	try(names(fit.se)<-paste(names(fit.final$par),".se",sep=""))
	#save(fit.final, file=paste(filename,'.fit.final',sep=""), compress=TRUE)
	#print(fit.final)
	#print(fit.se)
	final.matrix<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),fit.final$par),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu",names(fit.final$par))))
	#save(final.matrix, file=paste(filename,'.final.matrix',sep=""), compress=TRUE)
	rownames(final.matrix)<-paste("FINAL_",rownames(final.matrix),sep="") #to make it easier to grep
	print(formatC(final.matrix,format="f",digits=30,drop0trailing=TRUE))
	print(paste("FINAL_F ",F,sep=""))
	print(paste("FINAL_T ",T,sep=""))
	print(paste("FINAL_D ",D,sep=""))
	print(paste("FINAL_S ",S,sep=""))
	print(paste("FINAL_filename ",filename,sep=""))
	elapsedTime<-(proc.time()-startTime)[3]
	print(paste("elapsedTime = ",elapsedTime))
	final.matrix.ml<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),elapsedTime,coef(fit.final,full=TRUE,extra=TRUE)),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu","elapsedTime",names(coef(fit.final,full=TRUE,extra=TRUE)))))
	final.matrix.se<-matrix(c(NA,NA,NA,NA,NA,NA,NA,coef(fit.final.se,full=TRUE,extra=TRUE)),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu","elapsedTime",names(coef(fit.final.se,full=TRUE,extra=TRUE)))))
	final.matrix.all<-cbind(mle=final.matrix.ml,se=final.matrix.se)
	colnames(final.matrix.all)<-c("mle","se")
#	final.matrix.all<-matrix(c(fit.final$lnLik,AIC(fit.final,k=length(fit.final$par)),length(fit.final$par),length(grep("q",names(fit.final$par))),length(grep("lambda",names(fit.final$par))),length(grep("mu",names(fit.final$par))),coef(fit.final,full=TRUE,extra=TRUE),fit.se),ncol=1,dimnames=list(c("lnLik","AIC","k_all","k_q","k_lambda","k_mu",names(coef(fit.final,full=TRUE,extra=TRUE)),names(fit.se))))
	save(final.matrix.all, file=paste(filename,'.final.matrix.all',sep=""), compress=TRUE)
	#rownames(final.matrix.all)<-paste("FINALALL_",rownames(final.matrix.all),sep="") #to make it easier to grep
	#print(formatC(final.matrix.all,format="f",digits=30,drop0trailing=TRUE))


	#modify below this
	#diversificationVector<-getDiversificationRates(fit.final, type=diversificationType, partitionSize=S)
	#transitionVector<-getTransitionRates(fit.final,type=transitionType, partitionSize=S)
	#results.vector <- c(diversificationType, transitionType, logLik(fit.final), length(fit.final$par), getAIC(fit.final), diversificationVector, transitionVector)
	#save(results.vector,file=paste(filename,'raw.optim',sep=""),ascii=TRUE)
	#names(results.vector)<-c("diversificationType", "transitionType", "lnL", "K", "AIC", rep("div or trans ",length(results.vector)-5))
	#print(results.vector)
	#save(results.vector, file=paste(filename,'.optim',sep=""), ascii=TRUE)
	#shrink down final files
#	system("tail -40 run.Rout > tail.run.Rout")
#	system("/usr/bin/zip run.zip run.Rout")
#	system("rm run.Rout")
}


modify_transitions<-function(lik=lik, type=1, F=F, S=6, extralist=extralist) {
	#remember to see transition models in V*_UtilityFns.R
	maxStringLength=nchar(2^S) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.
	focalVector=F
	if (length(focalVector==1)) {
		focalVector<-stringToVector(focalVector)
	}



	#rather than typing manually all the restrictions, I will calculate this automatically
	constraintString="constrain(lik "
	for (charStateI in 1:((2^S))) {
		binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
		for (charStateJ in (charStateI+1):((2^S))) {
			if (charStateJ<=(2^S)) {
				binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
				numberMismatches=vectorMismatch(binaryStateIVector,binaryStateJVector) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
				if (numberMismatches>1) {
					constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),'~0, q',sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),'~0',sep="")
				}
			}
		}
	}

	#now for the actual model
	for (charStateI in 1:((2^S))) {
		binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
		fromFocal<-FALSE
		if( focalVector[comboAsDecimal(binaryStateIVector,S)] == 1) {
			fromFocal<-TRUE
		}
		for (charStateJ in 1:((2^S))) { #note we're starting from 1 here, too
			binaryStateJVector<-digitsBase(charStateJ-1,ndigits=S)[,1]
			numberMismatches=vectorMismatch(binaryStateIVector,binaryStateJVector) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
			if (numberMismatches==1) {
				toFocal<-FALSE
				if( focalVector[comboAsDecimal(binaryStateJVector,S)] == 1) {
					toFocal<-TRUE
				}
				qString<-"~qMost"
				if(fromFocal) {
					if (toFocal) {
						qString<-qFFbyModel[type] #see utility file
					}
					else {
						qString<-qFNbyModel[type]
					}
				}
				else {
					if (toFocal) {
						qString<-qNFbyModel[type]
					}
					else {
						qString<-qNNbyModel[type]
					}
				}
				constraintString=paste(constraintString,", q",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),sprintf(paste("%0",maxStringLength,"d",sep=""),charStateJ),qString,sep="")
			}
		}
	}

	replaceextralist(unlist(strsplit(extralistT[type],split=" ")))
	if (length(extralist)>0) {
		constraintString=paste(constraintString,", extra=c(",sep="")
		for (extraIndex in 1:length(extralist)) {
			constraintString=paste(constraintString,"'",extralist[extraIndex],"'",sep="")
			if (extraIndex<length(extralist)) {
				constraintString=paste(constraintString,", ",sep="")
			}
		}
		constraintString=paste(constraintString,")",sep="")
	}
	constraintString=paste(constraintString,")",sep="")
	print(paste("transition model: ",constraintString))
	return(eval(parse(text=constraintString)))
}

modify_diversification<-function(lik=lik, type=1, F=F, S=6, extralist=extralist) {
	print(paste("diversification input extralist = ",extralist))
	maxStringLength<-nchar(2^S) #assuming character states are single digits only works up to 2^3 states. If the max state is 64, diversitree counts 01, 02, etc.
	constraintString<-"constrain(lik "
	focalVector=F
	if (length(focalVector==1)) {
		focalVector<-stringToVector(focalVector)
	}

	for (charStateI in 1:((2^S))) {
		binaryStateIVector<-digitsBase(charStateI-1,ndigits=S)[,1]
		isFocal<-FALSE
		if( focalVector[comboAsDecimal(binaryStateIVector,S)] == 1) {
			isFocal<-TRUE
		}
		if (isFocal) {
			constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),bFbyModel[type],sep="")
			constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),dFbyModel[type],sep="")
		}
		else {
			constraintString=paste(constraintString,", lambda",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),bNbyModel[type],sep="")
			constraintString=paste(constraintString,", mu",sprintf(paste("%0",maxStringLength,"d",sep=""),charStateI),dNbyModel[type],sep="")
		}
	}
	origextralist<-extralist
	print(paste("extralistD is ",extralistD))
	replaceextralist(c(extralist, unlist(strsplit(extralistD[type],split=" "))))
	extralist<-c(origextralist, unlist(strsplit(extralistD[type],split=" ")))
	print(paste("extralist after replace extralist is ",extralist))
	if (length(extralist)>0) {
		constraintString=paste(constraintString,", extra=c(",sep="")
		for (extraIndex in 1:length(extralist)) {
			constraintString=paste(constraintString,"'",extralist[extraIndex],"'",sep="")
			if (extraIndex<length(extralist)) {
				constraintString=paste(constraintString,", ",sep="")
			}
		}
		constraintString=paste(constraintString,")",sep="")
	}
	constraintString=paste(constraintString,")",sep="")
	print(paste("diversification model: ",constraintString))
	return(eval(parse(text=constraintString)))
}

prepData<-function(P=P,F=F,T=T,D=D,S=6,sourcetraits="/data/abc/RunsJan2012/SourceData/Stebbins_prunenoper25i2012BCO.csv") {
	partitionVector<-strsplit(P,split="_")
	print(partitionVector)
	focalVector=F
	if (length(focalVector==1)) {
		focalVector<-stringToVector(focalVector)
	}

	partitionVector<-as.numeric(partitionVector[[1]])
	print(partitionVector)
	charsToInclude<-partitionVector[which(partitionVector>0)]
	stopifnot(S==length(charsToInclude))
	file<-sourcetraits
	phy<-"/data/abc/RunsJan2012/SourceData/floral_1.nex"
	tree<-read.nexus(phy)
	data<-read.csv(file)
	colnamesVector<-colnames(data)
	names<-colnamesVector[2:length(colnamesVector)]
	print(names)
	#keep only 1 (1st 25%) and 2 (2nd 25%) in Selection column
	#subdata<-subset(data,data$Selection<=2)
	#already done, so just
	subdata<-data #to minimize recoding

	#delete taxa with missing data anywhere
	for (colToExamine in 2:length(colnamesVector)) {
		subdata<-subdata[subdata[,colToExamine]!="?",]
	}

	#do recoding. Put chars at end of vector
	subdata[,(dim(subdata)[2]+1)]=subdata[,(charsToInclude[1]+1) ]
	if (length(charsToInclude)>=2) {
		for (charIndex in 2:length(charsToInclude)) {
			subdata[,(dim(subdata)[2])]=paste(subdata[,(dim(subdata)[2])], subdata[,(charsToInclude[charIndex]+1) ],sep="")
		}
	}
	#now make a new column that takes the 0, or 0 1, or 0 1 0 1, etc. and maps them into 1, 2, 3, 4
	#subdata[,(dim(subdata)[2]+1)]=1+todec(as.vector(unlist(strsplit(as.character(unlist(subdata[,(dim(subdata)[2]) ])," "))),mode="numeric"))
	subdata[,(dim(subdata)[2]+1)]<-NA
	#print(dim(subdata)[1])
	for (rowIndex in 1:dim(subdata)[1]) {
		rawData<-subdata[rowIndex,(dim(subdata)[2]-1)]
		rawData<-as.character(rawData)
		rawDataSplit<-strsplit(rawData,"")[[1]]
		rawDataSplit<-as.numeric(rawDataSplit)
		subdata[rowIndex,(dim(subdata)[2])]<-1+todec(rawDataSplit)
	}
	#if there are NA for the new char, this will take them out
	#if not, it will do nothing
	presentTaxa<-as.vector(subdata$Name_in_tree)
	to.drop <- setdiff(tree$tip.label, presentTaxa)
	tree2<-drop.tip(tree,to.drop)
	char<-subdata[,(dim(subdata)[2])]
	names(char)<-subdata$Name_in_tree
	colnamesVector<-colnames(data)
	print(colnamesVector)
	#finalname=colnamesVector[(charsToInclude[1]+1)]
	#print(finalname)
	#if (length(charsToInclude)>=2) {
	#	for (charIndex in 2:length(charsToInclude)) {
	#		finalname=paste(finalname, colnamesVector[(charsToInclude[charIndex]+1) ],sep="_")
	#	}
	#}
	finalname<-paste("T",T,"_D",D,"_F",F,sep="")

	#write out datafile for the character with matching tree
	write.csv(char,file=paste(finalname,".csv",sep=""))
	print(tree2)
	write.tree(tree2,file=paste(finalname,"tree",".t",sep=""))
	return(finalname)
}
