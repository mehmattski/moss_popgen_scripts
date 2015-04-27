# Given a file containing:
# Sporophytes (their IDs)
# Maternal Shoot ID (grouping)
# Microsatellite loci (pairs of columns)
# Size
# Mom and Dad Clone ID
# Percent Heterozygosity

# Calculate:

# FIS using a subsampled population approach (one sporophyte per maternal shoot)
# Effective number of alleles (Simpson's Index) and Shannon Diversity Index using subsampled population
# Paternity skew (Ke/K)
# Inbreeding depression (correlation between heterozygosity and size)

# The first two columns must be the individual IDs and the maternal shoot IDs.
# Then follows pairs of columns for the loci.
# After this, the column headers must be: size, momClone, dadClone, and het


#### READ IN THE FILE ###
setwd("examples")
filename="tenellum_2ns.txt"

rawdata = read.table(filename,header=T,sep="\t")

size = rawdata$size
momClone = rawdata$momClone
dadClone = rawdata$dadClone
het = rawdata$het

indata = subset(rawdata, select = -c(size,momClone,dadClone,het))

#### CREATE SUBSET POPULATIONS #######
numreps = 1000		#the number of subset populations


numMoms = length(unique(indata[,2]))	#Equals the number of samples in reduced datasets
numLoci = (length(indata) - 2)/2			#First two columns are IDs, rest is data

subsetrows = matrix(nrow=numMoms,ncol=numreps)	#Initialize subset matrix.
moms = levels((factor(indata[,2])))		#Get IDs of the moms

#Fill the subsetrows matrix with row index of one sample per mom, numreps times.
for(rep in 1:numreps){
	picks = c()
	for(x in moms){
		y = which(indata[,2]==x)
		pick = y[sample.int(length(y),1)] #Picks an index so when length(y) = 1, it picks that one sample.
		#pick = sample(y,1)
		picks = append(picks,pick)
	}
	subsetrows[,rep] = picks
}

indata = as.matrix(indata[,-c(1:2)])	#Matrix form of indata, lose the ID rows. Matrix goes faster than data frame. 

###### FIS, Ne, and I functions ########
obshet = function(x){
	#Input is two columns (one locus), output is observed heterozygosity.
	x[x==0]=NA
	x = na.omit(x)
	diffs = x[,1] - x[,2]
	numHets = nrow(x) - length(which(diffs==0)) 
	percentHet = numHets / nrow(x)
	return(percentHet)
	}
	
exphet = function(x){
	#Input is two columns (one locus), output is expected heterozygosity.
	#Calculated as the sum of the square of each allele's frequency.
	x[x==0]=NA
	x = x[!is.na(x)]	
	allele.counts = tabulate(x)
	he = 1 - sum((allele.counts/sum(allele.counts))^2)
	return (he)
	}	


get_fis = function(rep,loci){
	#Input is a vector of rows in original dataset, returns the mean fis for that rep.
	#loci is a vector of integers indicating which loci to be sampled.
	#in a bootstrap, loci might have repeated numbers
	
	rep.data = indata[rep,]	#Get rid of first two columns, matrix faster.

	#Initialize vectors for by-locus observed and expected heterozygosity and FIS
	obshet.values = rep(NA,numLoci)
	exphet.values = rep(NA,numLoci)
	for(locus in loci){
		column1 = locus*2 - 1		#Because each locus is two columns.
		column2 = locus*2			#Locus 3 = columns 5 and 6
		locusdata = rep.data[,column1:column2]				
		obshet.values[locus]=obshet(locusdata)	
		exphet.values[locus]=exphet(locusdata)
		}
	fis.values = 1 - obshet.values/exphet.values	#Calculate FIS for each subset from observed and expected heterozygosity.
	mean.fis = mean(fis.values,na.rm=T) #Remove values where exphet = 0
	return(mean.fis)
	}

allelefreqs = function(x){
	#Input is two columns (one locus), output is a vector of allele frequencies.
	table(x)/(nrow(x)*2)
}

simpson = function(rep,loci){
	rep.data = indata[rep,]	#Get rid of first two columns, matrix faster.

	#Initialize vectors for by-locus observed and expected heterozygosity and FIS
	ne = rep(NA,numLoci)
	for(locus in loci){
		column1 = locus*2 - 1		#Because each locus is two columns.
		column2 = locus*2			#Locus 3 = columns 5 and 6
		locusdata = rep.data[,column1:column2]				
		allele.freqs = allelefreqs(locusdata)
		ne[locus] = 1/sum(allele.freqs^2)	
		}
	mean.ne = mean(ne)
	return(mean.ne)
}

shannon = function(rep,loci){
	rep.data = indata[rep,]	#Get rid of first two columns, matrix faster.

	#Initialize vectors for by-locus observed and expected heterozygosity and FIS
	i = rep(NA,numLoci)
	for(locus in loci){
		column1 = locus*2 - 1		#Because each locus is two columns.
		column2 = locus*2			#Locus 3 = columns 5 and 6
		locusdata = rep.data[,column1:column2]				
		allele.freqs = allelefreqs(locusdata)
		i[locus] = -1*sum(allele.freqs*log(allele.freqs))	
		}
	mean.i = mean(i)
	return(mean.i)
}
###### CALCULATE STUFF #########

fis.values = apply(subsetrows,2,get_fis,c(1:numLoci))
fis.boot95 = quantile(fis.values,c(0.025,0.975))
effective.alleles = apply(subsetrows,2,simpson,c(1:numLoci))
simpson.boot95 = quantile(effective.alleles,c(0.025,0.975))
shannon.indicies = apply(subsetrows,2,shannon,c(1:numLoci))
shannon.boot95 = quantile(shannon.indicies,c(0.025,0.975))

##### Inbreeding Depression #####
inbreeding.lm = lm(size~het)
inbreeding.summary = summary(inbreeding.lm)
inbreeding.slope = inbreeding.lm$coefficients[2]
inbreeding.rsquared = inbreeding.summary$r.squared
inbreeding.pvalue = inbreeding.summary$coefficients[2,4]


#### Multiple Paternity and Paternity Skew ######
dads.by.momshoot = split(dadClone,rawdata$momID)
k.by.momshoot = unlist(lapply(dads.by.momshoot,function(x) length(unique(x))))
prop.dad = lapply(dads.by.momshoot, function(x) table(x)/length(x))
ke.by.momshoot = unlist(lapply(prop.dad,function(x) 1/sum(x^2))) 
more.than.one.dad = which(k.by.momshoot > 1)
# For Ke/K, remove shoots that were only fertlized by one dad
k.singletons.removed = k.by.momshoot[more.than.one.dad]
ke.singletons.removed = ke.by.momshoot[more.than.one.dad]
ke.k.by.momshoot= ke.singletons.removed/k.singletons.removed

k = mean(k.by.momshoot)
ke = mean(ke.by.momshoot)
ke.k = mean(ke.singletons.removed/k.singletons.removed)
# Bootstrap K, Ke, and Ke/K to see if > 1
numBoots = 100
k.resamples = lapply(1:numBoots,function(i) sample(k.by.momshoot,replace=T))
k.boot = unlist(lapply(k.resamples,mean))
ke.resamples = lapply(1:numBoots,function(i) sample(ke.by.momshoot,replace=T))
ke.boot = unlist(lapply(ke.resamples,mean))
ke.k.resamples = lapply(1:numBoots, function(i) sample(ke.k.by.momshoot,replace=T))
ke.k.boot = unlist(lapply(ke.k.resamples,mean))






####Put it together ####
all_stats = c(
	nrow(rawdata),
	length(unique(rawdata$momID)),
	mean(effective.alleles),
	simpson.boot95,
	mean(shannon.indicies),
	shannon.boot95,
	mean(fis.values),
	fis.boot95,
	inbreeding.slope,
	inbreeding.rsquared,
	inbreeding.pvalue,
	k,
	quantile(k.boot,c(0.025,0.975)),
	ke,
	quantile(ke.boot,c(0.025,0.975)),
	ke.k,
	quantile(ke.k.boot,c(0.025,0.975))
)

#For each stat, .lo and .hi refer to the 2.5 and 97.5 percentile within each pseudopopulation.

names(all_stats) = c(
	"num2ns", #Number of Sporophytes
	"numMoms", #Number of Mothers
	"Ne",	#Effective Number of Alleles
	"Ne.lo", 
	"Ne.hi", 
	"I",	#Shannon's Allelic Diversity 
	"I.lo", 
	"I.hi",
	"FIS", #Inbreeding coefficient
	"FIS.lo",
	"FIS.hi",
	"ID.slope", #Inbreeding depression slope
	"ID.r2",	#Inbreeding depression r-squared
	"ID.p",		#P-value for Inbreeding depression slope
	"k",		#Multiple paternity
	"k.lo",
	"k.hi",
	"ke",		#Effective number of parents
	"ke.lo",
	"ke.hi",
	"ke.k",		#Paternity Skew
	"ke.k.lo",
	"ke.k.hi"
)

#Write all_stats values to the clipboard
clip = pipe("pbcopy",'w')
write.table(all_stats,file=clip,row.names=F,col.names=F)
close(clip)

