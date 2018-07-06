# ----------------------------------------------------------------------------------
# David Phillips
#
# 7/3/2018
# Prep covariates and their relationship with vaccine coverage at the pixel level
# ----------------------------------------------------------------------------------


# --------------------
# Set up R
rm(list=ls())
library(boot)
library(stats)
library(data.table)
library(plyr)
library(raster)
library(rgeos)
library(lme4)
# --------------------


# ---------------------------------------------------
# Parameters and settings

# whether to decrease the pixel resolution for speed 
decrease_res = TRUE
# ---------------------------------------------------


# --------------------------------------------------------------------------------
# Files and directories

j = ifelse(Sys.info()[1]=='Windows', 'J:', '/home/j')

# directory for covariates
covDir = paste0(j, '/WORK/11_geospatial/01_covariates/00_MBG_STANDARD/')

# covariate files
eduFile = paste0(covDir, 'edu_mean_stage_2/mean/1y/edu_mean_stage_2_mean_1y_2016_00_00.tif')
accessFile = paste0(covDir, 'access2/mean/synoptic/access2_mean_synoptic.tif')

# directory for vaccine estimates
vacDir = '/share/geospatial/mbg/vaccine/dpt3_cov/output/2018_06_11_22_17_24/'

# vaccine coverage file
vacFile = paste0(vacDir, 'dpt3_cov_mean_raked_raster.tif')

# location for admin1 shapefile
admin1Dir = paste0(j, '/WORK/11_geospatial/06_original shapefiles/GADM_3_6/final2/')

# shapefile files
admin0File = paste0(admin1Dir, 'gadm36_0.shp')
admin1File = paste0(admin1Dir, 'gadm36_1.shp')

# directory for output
outDir = paste0(j, '/temp/davidp6/vaccine_coverage_correlates/')

# output files
outFile = paste0(outDir, 'input_data.rdata')
# --------------------------------------------------------------------------------


# -----------------------------------------------------
# Load/prep data

# load covariates
edu = raster(eduFile)
access = raster(accessFile)

# load estimates
vac = raster(vacFile)

# project covariates to estimate pixels to ensure merge
edu = projectRaster(edu, vac)
access = projectRaster(access, vac)

# crop to remove white space
edu = crop(edu, extent(vac))
access = crop(access, extent(vac))

# mask out countries with no coverage estimates
edu = mask(edu, vac)
access = mask(access, vac)

# aggregate resolution if specified
if (decrease_res) { 
	edu = aggregate(edu, fact=8)
	access = aggregate(access, fact=8)
	vac = aggregate(vac, fact=8)
}

# convert to data tables
eduDT = data.table(as.data.frame(edu, xy=TRUE))
accessDT = data.table(as.data.frame(access, xy=TRUE))
vacDT = data.table(as.data.frame(vac, xy=TRUE)) 

# merge pixels
data = merge(eduDT, accessDT, by=c('x','y'),all=TRUE)
data = merge(data, vacDT, by=c('x','y'),all=TRUE)

# make cell id to identify polygons later
cellNumbers = cellFromXY(vac, data[,c('x','y'), with=FALSE])
data[, cell:=cellNumbers]

# omit na to reduce extent
data = na.omit(data)

# load/simplify shapefile at country level
map0 = shapefile(admin0File)
map0tmp = map0@data
# map0 = gSimplify(map0, tol=.1, topologyPreserve=TRUE)
# map0 = as(map0, 'SpatialPolygonsDataFrame')
# map0@data = map0tmp
map0 = crop(map0, extent(vac))

# load/simplify shapefile at admin-1 level
map1 = shapefile(admin1File)
map1tmp = map1@data
# map1 = gSimplify(map1, tol=.1, topologyPreserve=TRUE)
# map1 = as(map1, 'SpatialPolygonsDataFrame')
# map1@data = map1tmp
map1 = crop(map1, extent(vac))
# -----------------------------------------------------


# ----------------------------------------------------------------------------------
# Identify the admin1 for each pixel #

# extract pixels for each polygon
extractionOrig = extract(vac, map1, cellnumbers=TRUE) # SLOW
extraction = extractionOrig

# identify the polygon of each pixel from the extraction
for(i in seq(length(extraction))) { 
	if (!is.null(ncol(extraction[[i]]))) colnames(extraction[[i]]) = c('cell','value')
	extraction[[i]] = cbind(extraction[[i]], i)
}
extraction = data.table(do.call('rbind.fill.matrix',extraction))

# convert ids to names
names = data.table(map1@data[,c('NAME_1','GID_1')]) 
setnames(names, c('admin1_name','admin1_id'))
names[, i:=seq(.N)]
extraction = merge(extraction, names, 'i')

# store ids and names alone
ids = extraction[, c('cell','admin1_name', 'admin1_id'), with=FALSE]

# identify which boundary each pixel is in
data = merge(data, ids, 'cell')
# ----------------------------------------------------------------------------------


# ------------------------------------------------------
# Simplify maps
map0tmp = map0@data
map0 = gSimplify(map0, tol=.1, topologyPreserve=TRUE)
map0 = as(map0, 'SpatialPolygonsDataFrame')
map0@data = map0tmp
map1tmp = map1@data
map1S = gSimplify(map1, tol=.1, topologyPreserve=TRUE)
map1S = as(map1S, 'SpatialPolygonsDataFrame')
map1S@data = map1tmp
# ------------------------------------------------------


# --------------------------------------------------------------------------------------------------------
# Run analyses

# estimate overall explained variance with a linear model
lmFit = lm(logit(dpt3_cov_mean_raked_raster) ~ log(access2_mean_synoptic) + edu_mean_stage_2_mean_1y_2016_00_00, data=data) 

# print explained variance
af = anova(lmFit)
afss = af$'Sum Sq'
explainedVariancesOverall = data.table(variable=c('access','edu','supply'), explained_variance=afss/sum(afss))
explainedVariancesOverall

# estimate admin1-level random effects
lmeFit = lmer(logit(dpt3_cov_mean_raked_raster) ~ log(access2_mean_synoptic) + edu_mean_stage_2_mean_1y_2016_00_00 + 
				(1 + edu_mean_stage_2_mean_1y_2016_00_00 | admin1_id) + (1 + (access2_mean_synoptic) | admin1_id), data=data) 

# extract admin1-level variances
V1<-VarCorr(lmeFit)$admin1_id
V2<-VarCorr(lmeFit)$admin1_id.1

# compute variance explained at each admin1 mean
means = data[,.(mean(edu_mean_stage_2_mean_1y_2016_00_00), mean(log(access2_mean_synoptic))), by='admin1_id']
Z1<-cbind(rep(1,nrow(means)), means[[2]])
Z2<-cbind(rep(1,nrow(means)), means[[3]])
explainedVariance1<-diag(Z1%*%V1%*%t(Z1))
explainedVariance2<-diag(Z2%*%V2%*%t(Z2))

# convert explained variances to proportions
resid = sqrt(attr(VarCorr(lmeFit), 'sc')) # IS THERE A WAY TO GET GROUP-SPECIFIC RESIDUAL VARIANCE ESTIMATES?
explainedVarianceSc1 = explainedVariance1/(explainedVariance1+explainedVariance2+resid)
explainedVarianceSc2 = explainedVariance2/(explainedVariance1+explainedVariance2+resid)
explainedVariancesAdmin1 = data.table(admin1_id=accessMeans$admin1_id, edu=explainedVarianceSc1, access=explainedVarianceSc2)

# store admin1-level explained variance for each factor
explainedVariancesAdmin1[, supply:=1-(edu+access)]
explainedVariancesAdmin1
explainedVariancesAdmin1[,lapply(.SD, mean), .SDcols=c('edu','access','supply')]
explainedVariancesAdmin1[,lapply(.SD, sd), .SDcols=c('edu','access','supply')]
# --------------------------------------------------------------------------------------------------------


# --------------------------------------------------------
# Save
objs = NULL
for (o in c('access','edu','vac','data','map0','map1','lmFit','lmeFit','explainedVariancesOverall','explainedVariancesAdmin1') {
	if (o %in% ls()) objs = c(objs, o)
}
save(list=objs, file=outFile)
# -------------------------------------------------------
