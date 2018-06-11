library(diversitree) #obvious
library(sfsmisc) #for counting in binary
library(partitions) #for converting from binary back to decimal
library(gmp) #for dealing with big integers


#' basic idea:
#' 6 binary traits
#' means 2^6 = 64 possible character-state combinations (note: not all of them are realized)
#' select between 0 and 32 of these combos as "focal" ones
#' syntax for focal set: a 64-long binary string, telling whether a given combo is in the focal set
#' filter these for meaningful sets: a set consisting of combos that are 0*1**** is meaningful, for example
#' [problem: equal sets: 0****** vs 1****** are the same under full models (i.e., different birth, death, and transition rates -- but only in that case]

#' transition models:
#'	There are four transition rates q_focal->nonfocal, q_focal->focal, q_nonfocal->nonfocal, q_nonfocal->focal = qFN, qFF, qNF, qNF. Note that you can have within set transitions where the number of combos in a set is greater than 1: set 0*1**** has a 0010000->0010001 transition rate
#'	1: all rates the same: qNF=qNN=qFF=qFN
#'	2: inflow different: qNF vs qNN=qFF=qFN
#'	3: outflow different: qFN vs qNN=qFF=qNF
#'	4: inflow and outflow different: qFN vs qNF vs qFF=qNN
#'	5: freedom: qFN vs qNF vs qFF vs qNN

#' Get the transition models
#'
#' Simple function to generate the potential transition models. Modify as you wish
#'
#' @return data.frame of transition models
#' @export
transitionModels <- function() {
	transitionModels<-data.frame(cbind(c(1:5),c(1,2,2,3,4),c(0,1,1,1,2),c("equal","inflow","outflow","inandoutflow","free")),stringsAsFactors=FALSE)
	names(transitionModels)<-c("T","k_q","min_focalcombos","description") #the min focal states is because some models don't make sense in certain cases. For example, if you have one focal state, models that have qFF as a free parameter don't make sense
	return(transitionModels)
}

qFFbyModel <- function() {
	return(c("~qMost", "~qMost", "~qMost", "~qMost", "~qFF"))
}


qNNbyModel <- function() {
	return(c("~qMost", "~qMost", "~qMost", "~qMost", "~qNN"))
}


qNFbyModel <- function() {
	return(c("~qMost", "~qNF"  , "~qMost", "~qNF"  , "~qNF"))
}

qFNbyModel <- function() {
	return(c("~qMost", "~qMost", "~qFN"  , "~qFN"  , "~qFN"))
}


extralistT <- function() {
	return(c("qMost", "qMost qNF", "qMost qFN", "qMost qNF qFN", "qFF qNN qNF qFN")) #stores the extralist values; remember to do strsplit before use
}


# note that if focal set size is 1, qFF is not a free parameter, so change K accordingly

# diversification models:
#	There are four diversification rates: bF, bN, dF, dN
#	1: yule: bF=bN, dF=dN=0
#	2: two-rate yule: bF vs bN vs dF=dN=0
#	3: simple birth-death: bF=bN, dF=dN
#	4: two birth, one death: bF vs bN vs dF=dN
#	5: one birth, two death: bF=bN vs dF vs dN
#	6: freedom: bF vs bN vs dF vs dN

diversificationModels <- function() {
	data.frame(cbind(c(1:6),c(1,2,1,2,1,2),c(0,0,1,1,2,2),c(0,1,0,1,1,1),c("yule","tworateyule","birthdeath","twobirthonedeath","onebirthtwodeath","free")),stringsAsFactors=FALSE)
	names(diversificationModels)<-c("D","k_b","k_d","min_focalcombos","description")
	return(diversificationModels)
}

bFbyModel <- function() {
	return(c("~lambdaAll", "~lambdaF",  "~lambdaAll", "~lambdaF", "~lambdaAll", "~lambdaF"))
}


bNbyModel <- function() {
	return(c("~lambdaAll", "~lambdaN",  "~lambdaAll", "~lambdaN", "~lambdaAll", "~lambdaN"))
}

dFbyModel <- function() {
	return(c("~0"        , "~0"      ,  "~muAll"    , "~muAll"  , "~muF"      , "~muF"    ))
}

dNbyModel <- function() {
	return(c("~0"        , "~0"      ,  "~muAll"    , "~muAll"  , "~muN"      , "~muN"    ))
}

extralistD <- function() {
	return(c("lambdaAll", "lambdaF lambdaN", "lambdaAll muAll", "lambdaF lambdaN muAll", "lambdaAll muF muN", "lambdaF lambdaN muF muN"))
}

# Definitions:
#   character: single trait, like petal symmetry
#	character combination: combination of traits, like 0100011
#   focal set: group of character combinations, like c(0100011, 0100010) is one possible focal set (in this case 010001*)

# max number of parameters in this case is just 8, not bad for 500 taxa. Want to only do one case where all states are equal for all rates

# constrain root state = 0000000

#Global definitions
nchar=6
S=nchar
partitionSize<-nchar

#utility functions
toBinLarge<-function (x, base = 2, S=6)  #modification from sfsmisc package to deal with large numbers
{
	ndigits = 2^S
	if (class(x)!="bigz") {
		x<-as.bigz(x)
	}
	if (class(base)!="bigz") {
		base<-as.bigz(base)
	}
	binaryVector<-c()
	while (x>0) {
		binaryVector<-c(as.numeric(x%%base),binaryVector)
		x<-x%/%base
	}
	fullVector<-rep(0,ndigits)
	lengthDiff<-length(fullVector)-length(binaryVector)
	if (lengthDiff<0) {
		lengthDiff<-0
		fullVector<-rep(0,length(binaryVector))
	}
	if (length(binaryVector)>0) {
		for (i in 1:length(binaryVector)) {
			fullVector[i+lengthDiff]<-binaryVector[i]
		}
	}
	return(fullVector)
}


vectorMismatch<-function(vector1, vector2) {
	if (length(vector1)!=length(vector2)) {
		return(NA)
	}
	else {
		return(sum(1-(vector1==vector2))) #so we have two vectors, say 00101 and 00110 (though as length 5 vectors). Doing v1==v2 leads to T T T F F. T=1 for R and F=0, so 1-(v1==v2) = c(1-1,1-1,1-1,1-0,1-0), sum of which is the number of mismatches
	}
}

vectorMismatchExcludePositions<-function(vector1, vector2, excludePositions) {
	if (length(vector1)!=length(vector2)) {
		return(NA)
	}
	else if (length(excludePositions)==0) {
		return(vectorMismatch(vector1, vector2))
	}
	else if (max(excludePositions)>max(length(vector1),length(vector2))) {
		return(NA)
	}
	else {
		return(sum(1-((vector1==vector2)[-1*excludePositions])))
	}
}

#comboDecimal go from 1:2^S
comboAsBinaryVector<-function(comboDecimal,S=6) { #works best of combo Decimal is bigz class
	if (class(comboDecimal)!="bigz") {
		comboDecimal<-as.bigz(comboDecimal)
	}
	return(toBinLarge(comboDecimal-1,S=S)) #since this goes from 000000 -> 111111
}

comboAsBinaryString<-function(comboDecimal,S=6) {
	return(paste(comboAsBinaryVector(comboDecimal,S),sep="",collapse=""))
}

comboAsDecimal<-function(comboBinary,S=6) {
	if (length(comboBinary)==1) { #we have a string, first make it a vector
		comboBinary<-strsplit(comboBinary,split="")[[1]]
	}
	return(1+todec(as.numeric(comboBinary)))
}



#focalAsBinaryVector goes from 000......000 (nothing is special) to 111......111 (all is special). Is of length 2^6.
#so if it is 110000...000  only combo numbers 1 and 2 are in the focal set
focalAsBinaryVector<-function(focalDecimal,S=6) { #works best if focalDecimal is bigz class
	if (class(focalDecimal)!="bigz") {
		warning("converting from a decimal to a focal binary vector will not work well if the decimal is not in bigz format for large values (R
with integers)")
	}
	return(toBinLarge(focalDecimal,S=S))
}

focalAsBinaryString<-function(focalDecimal,S=6) {
	return(paste(focalAsBinaryVector(focalDecimal,S),sep="",collapse=""))
}

focalAsDecimal<-function(focalBinary,S=6) {
	if (length(focalBinary)==1) { #we have a string, first make it a vector
		focalBinary<-strsplit(focalBinary,split="")[[1]]
	}
	return(todec(as.numeric(focalBinary)))
}

convertFocalToCombos<-function(focalBinaryVector) {
	return(which(focalBinaryVector==1))
}

numberFocalCombos<-function(focalBinaryVector) {
	return(length(which(focalBinaryVector==1)))
}


stringToVector<-function(inString) {
	return(as.numeric(unlist(strsplit(inString,split=""))))
}

vectorToString<-function(inVector) {
	return(paste(inVector,sep="",collapse=""))
}

maxFocalAsBinaryVector<-rep(1,2^S)


createComboMatrix<-function(focalBinaryVector,S=6) {
	focalCombos<-convertFocalToCombos(focalBinaryVector)
	comboMatrix<-matrix()
	if (length(focalCombos)>0) {
		for (i in 1:length(focalCombos)) {
			if (i==1) {
				comboMatrix<-matrix(comboAsBinaryVector(as.bigz(focalCombos[i]),S),nrow=1)
			}
			else {
				comboMatrix<-rbind(comboMatrix,comboAsBinaryVector(as.bigz(focalCombos[i]),S))
			}
		}
	}
	return(comboMatrix)
}

getFocalSummaryLabel<-function(focalBinaryVector,S=6,any="*") {
	comboMatrix<-createComboMatrix(focalBinaryVector,S)
	labelVector<-rep(any,S)
	if (numberFocalCombos(focalBinaryVector)>1) {
		for (i in 1:dim(comboMatrix)[2]) {
			if (var(comboMatrix[,i])==0) {
				labelVector[i]<-as.character(comboMatrix[1,i])
			}
		}
	}
	else if (numberFocalCombos(focalBinaryVector)==1) {
		labelVector<-as.character(comboMatrix[1,])
	}
	return(labelVector)
}

interestingFocal<-function(focalBinaryVector,S=6) {
	interestingFocal<-FALSE
	if (sum(focalBinaryVector)>0.5*length(focalBinaryVector)) { #too many focal combos (don't want more than half the combos)
		return(FALSE)
	}
	comboMatrix<-createComboMatrix(focalBinaryVector,S)
	focalCombos<-convertFocalToCombos(focalBinaryVector)
	if (length(focalCombos)<2) {
		interestingFocal<-TRUE #is interesting because a single focal trait or zero focal traits
	}
	else {
		numberInvariantSites<-length(which(colMeans(comboMatrix)==0)) + length(which(colMeans(comboMatrix)==1)) #these are the invariant columns
		if ((2^(S-numberInvariantSites))==dim(comboMatrix)[1]) {
			interestingFocal<-TRUE #imagine just 3 chars. c(010,011) works as 01*, but c(010,111) does not (it is just *1*, but only part of *1*, omitting 110 and 011). This tests that.
		}
	}
	return(interestingFocal)
}

getAllInterestingFocalVectorsInefficient<-function(S=6, file="interestingFocalDataFrame.Rsave"){
	focalDecimal<-as.bigz(0)
	focalVector<-focalAsBinaryVector(focalDecimal,S)
	totalInterestingFocal<-0
	totalAllFocal<-0
	interestingFocalDataFrame<-data.frame()
	while (min(focalVector)==0) { #so this will stop once focalVector gets to 111...1111
		totalAllFocal<-totalAllFocal+1
		if (interestingFocal(focalVector,S)) { #if this is an interesting combination
			totalInterestingFocal<-totalInterestingFocal+1
			if (totalInterestingFocal==1) {
				interestingFocalDataFrame<-data.frame(cbind(as.character.bigz(focalDecimal),vectorToString(getFocalSummaryLabel(focalVector,S,"x")),numberFocalCombos(focalVector),vectorToString(focalVector)),stringsAsFactors=FALSE)
			}
			else {
				interestingFocalDataFrame<-rbind(interestingFocalDataFrame,data.frame(cbind(as.character.bigz(focalDecimal),vectorToString(getFocalSummaryLabel(focalVector,S,"x")),numberFocalCombos(focalVector),vectorToString(focalVector)),stringsAsFactors=FALSE))
			}
			save(interestingFocalDataFrame,file=file,compress=TRUE)
			print(paste(vectorToString(getFocalSummaryLabel(focalVector,S)), ":",totalInterestingFocal, "/",totalAllFocal,vectorToString(focalVector),sep=" ",collapse=""))
		}
		focalDecimal<-focalDecimal+1
		focalVector<-focalAsBinaryVector(focalDecimal,S)
	}
}

convertFocalLabelToFocalVector<-function(focalLabel,S=6,uncertainty="2") {
	labelVector<-strsplit(focalLabel,split="")[[1]]
	uncertainChars<-which(labelVector==uncertainty)
	focalVector<-rep(0,2^S)
	for (i in 1:2^S) {
		if(0==vectorMismatchExcludePositions(labelVector,comboAsBinaryVector(as.bigz(i),S),uncertainChars)) {
			focalVector[i]<-1
		}
	}
	return(focalVector)
}


#Here's the new one that includes 2******
getAllInterestingFocalVectorsStringsEfficient<-function(S=6) {
  focalVectorList<-list()
   maxNumber<-3^S #could have state 0, 1, or 2
   for (i in 0:(maxNumber-1)) {
      #print(paste("now working on branch ",i))
      focalVectorLabel<-vectorToString(toBinLarge(i,3,S))
     # print(focalVectorLabel)
      focalVector<-convertFocalLabelToFocalVector(focalVectorLabel,S,"2")
      #print(focalVector)
#      names(focalVector)<-focalVectorLabel
#      print(focalVectorLabel)
      focalVectorList<-append(focalVectorList,vectorToString(focalVector))
      names(focalVectorList[[i+1]])<-focalVectorLabel
     # print(focalVectorList)
   }
   return(focalVectorList)
}

summarizeModelWeights<-function(summary.dataframe=summary.dataframe,S=6,transitionModels=transitionModels, diversificationModels=diversificationModels) {
	modelWeights<-matrix(0,nrow=1+dim(transitionModels)[1],ncol=1+dim(diversificationModels)[1])
	for (transitionIndex in 1:dim(transitionModels)[1]) {
		t.summary.dataframe<-subset(summary.dataframe,T==transitionIndex)
		modelWeights[transitionIndex,1+dim(diversificationModels)[1] ] <-sum(t.summary.dataframe$AICweight)
		for (diversificationIndex in 1:dim(diversificationModels)[1]) {
			td.summary.dataframe<-subset(t.summary.dataframe,D==diversificationIndex)
			modelWeights[transitionIndex,diversificationIndex]<-sum(td.summary.dataframe$AICweight)
			modelWeights[1+dim(transitionModels)[1],diversificationIndex]<- modelWeights[1+dim(transitionModels)[1],diversificationIndex] + modelWeights[transitionIndex,diversificationIndex]
		}
	}
	modelWeights[1+dim(transitionModels)[1],1+dim(diversificationModels)[1] ]<-1
	dimnames(modelWeights)<-list(c(transitionModels$description,"diversification totals"),c(diversificationModels$description,"transition totals"))
	print(round(modelWeights),digits=3)
	return(modelWeights)
}


#this is based on the implementation in lmomco by William H. Asquith
#His code incorporates the correction used by EPA for zero values
#I have just added weights to it
weightedHarmonicMeanZeroCorrection <- function(x, w=rep(1,length(x))) {
	w <- w[!is.na(x)]
	w <- w/sum(w)
	x <- x[!is.na(x)]
    n <- length(x)
    x.nonzero <- x[x != 0]
    w.nonzero <- w[x != 0]
    w.zero <- w[x == 0]
    n.zero <- n - length(x.nonzero)
    HM <- weighted.mean(c(1/weighted.mean(1/x.nonzero, w.nonzero), 0), c(sum(w.nonzero), sum(w.zero)))
    if(n.zero == n) {
    	HM <- 0
    }
    return(HM)

}



getD.ef<-function(ef, turn) {
  return(ef*turn/(1+ef))
}

getB.ef<-function(ef, turn) {
  return(turn-getD(ef, turn))
}

getD.net<-function(net.div, turn) {
  return((turn - net.div)/2)
}

getB.net<-function(net.div, turn) {
  return(net.div + getD.net(net.div, turn))
}
