# ----------------------------------------------------------------------------------
# David Phillips
#
# 7/11/2018
# Analyze covariates and their relationship with vaccine coverage at the pixel level
# Intended to be run by prep_correlates.r. 
# ( It really just moves some more complex code out of that script)
# ----------------------------------------------------------------------------------


# --------------------
# Set up R
library(boot)
library(stats)
library(data.table)
library(plyr)
library(raster)
library(rgeos)
library(stringr)
library(lme4)
# --------------------


# --------------------------------------------------------------------------------------------------------
# Run analyses

# create additional variables for analysis
pctle = quantile(data[access2_mean_synoptic!=0]$access2_mean_synoptic, .025)
data[, access_offset:=access2_mean_synoptic]
data[, access_offset:=access_offset+pctle]
data[, iso3:=str_sub(admin1_id, 1, 3)]

# estimate overall explained variance with a linear model
lmFit = lm(logit(dpt3_cov_mean_raked_raster) ~ log(access_offset) + edu_mean_stage_2_mean_1y_2016_00_00, data=data) 

# estimate admin1-level random effects
lmeFit = lmer(logit(dpt3_cov_mean_raked_raster) ~ log(access_offset) + edu_mean_stage_2_mean_1y_2016_00_00 + 
				(1 + edu_mean_stage_2_mean_1y_2016_00_00 + log(access_offset) | admin1_id), data=data)
# --------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------
# Assess explained variance

# print explained variance from simple model
af = anova(lmFit)
afss = af$'Sum Sq'
explainedVariancesOverall = data.table(variable=c('access','edu','supply'), explained_variance=afss/sum(afss))
explainedVariancesOverall


# extract admin1-level variances if independence assumed between random effects (not the perfect way to branch)
if (!is.null(VarCorr(lmeFit)$admin1_id.1)) {  
	V1<-VarCorr(lmeFit)$admin1_id
	V2<-VarCorr(lmeFit)$admin1_id.1

	# compute variance explained at each admin1 mean
	# ref: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q3/002688.html
	means = data[,.(mean(edu_mean_stage_2_mean_1y_2016_00_00), mean(log(access_offset))), by='admin1_id']
	Z1<-cbind(rep(1,nrow(means)), means[[2]])
	Z2<-cbind(rep(1,nrow(means)), means[[3]])
	explainedVariance1<-diag(Z1%*%V1%*%t(Z1))
	explainedVariance2<-diag(Z2%*%V2%*%t(Z2))

	# compute proportion of explained variance
	resid = sqrt(attr(VarCorr(lmeFit), 'sc'))
	explainedVarianceSc1 = explainedVariance1/(explainedVariance1+explainedVariance2+resid)
	explainedVarianceSc2 = explainedVariance2/(explainedVariance1+explainedVariance2+resid)
	explainedVariancesAdmin1 = data.table(admin1_id=means$admin1_id, edu=explainedVarianceSc1, access=explainedVarianceSc2)

	# store admin1-level explained variance for each factor
	explainedVariancesAdmin1[, supply:=1-(edu+access)]
	explainedVariancesAdmin1
	explainedVariancesAdmin1[,lapply(.SD, mean), .SDcols=c('edu','access','supply')]
	explainedVariancesAdmin1[,lapply(.SD, sd), .SDcols=c('edu','access','supply')]
	
# extract admin1-level variances if correlation between random effects	
} else { 

	# extract admin1-level variances
	V1<-VarCorr(lmeFit)$admin1_id

	# set up design matrix at the mean of each admin1
	means = data[,.(edu_mean_stage_2_mean_1y_2016_00_00=mean(edu_mean_stage_2_mean_1y_2016_00_00), access_offset=mean(log(access_offset))), by='admin1_id']
	Z1<-cbind(rep(1,nrow(means)), means$edu_mean_stage_2_mean_1y_2016_00_00, means$access_offset)
	Z2<-cbind(rep(1,nrow(means)), means$access_offset, means$edu_mean_stage_2_mean_1y_2016_00_00)

	# compute variance explained at each admin1 mean
	# ref: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q3/002688.html
	explVar1<-diag(Z1%*%V1%*%t(Z1))
	explVar2<-diag(Z2%*%V1%*%t(Z2))

	# extract unexplained variance overall
	residVar = sqrt(attr(VarCorr(lmeFit), 'sc'))

	# compute proportion of explained variance as a simple fraction
	explVarPct1 = explVar1/(explVar1+explVar2+residVar)
	explVarPct2 = explVar2/(explVar1+explVar2+residVar)
	explainedVariancesAdmin1 = data.table(admin1_id=means$admin1_id, edu=explVarPct1, access=explVarPct2)
	
	# store admin1-level explained variance for each factor
	explainedVariancesAdmin1[, supply:=1-(edu+access)]
	explainedVariancesAdmin1
	explainedVariancesAdmin1[,lapply(.SD, mean), .SDcols=c('edu','access','supply')]
	explainedVariancesAdmin1[grepl('UGA',admin1_id),lapply(.SD, mean), .SDcols=c('edu','access','supply')]
	explainedVariancesAdmin1[grepl('ZMB',admin1_id),lapply(.SD, mean), .SDcols=c('edu','access','supply')]
	explainedVariancesAdmin1[,lapply(.SD, sd), .SDcols=c('edu','access','supply')]
}
# --------------------------------------------------------------------------------------------------------
