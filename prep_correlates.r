# ----------------------------------------------------------------------------------
# David Phillips
#
# 7/3/2018
# Prep covariates and their relationship with vaccine coverage at the pixel level
# The current working directory should be the root of this repo
# Inputs:
# decrease_res - logical. whether to decrease the resolution of the shapefiles and rasters for size/speed
# res_factor - numeric. a factor by which to decrease raster resolution
# subset - character. either a vector of ISO3 codes, "Africa" or "All"
# Outputts:
# nothing. saves numerous objects to the path outFile, defined by the input parameters
# ----------------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Define function
prepCorrelates = function(decrease_res=TRUE, res_factor=8, subset='Africa') { 
# ----------------------------------------------------------------------------


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

	# shapefile of lakes
	shapeFileLakes = paste0(j, '/WORK/11_geospatial/06_original shapefiles/GLWD_lakes/glwd_1.shp')

	# directory for output
	outDir = paste0(j, '/temp/davidp6/vaccine_coverage_correlates/data')

	# output files
	outFile = paste0(outDir, 'input_data_', paste0(subset, collapse=''), 'resolution', res_factor, '.rdata')

	# preset for all of africa for convenience
	if (subset[1]=='Africa') subset = c('DZA', 'AGO', 'SHN', 'BEN', 'BWA', 'BFA', 'BDI', 'CMR', 'CPV', 'CAF', 'TCD', 'COM', 'COG', 'COD', 'DJI', 'EGY', 'GNQ', 'ERI', 'SWZ', 'ETH', 'GAB', 'GMB', 'GHA', 'GIN', 'GNB', 'CIV', 'KEN', 'LSO', 'LBR', 'LBY', 'MDG', 'MWI', 'MLI', 'MRT', 'MUS', 'MYT', 'MAR', 'MOZ', 'NAM', 'NER', 'NGA', 'STP', 'REU', 'RWA', 'STP', 'SEN', 'SYC', 'SLE', 'SOM', 'ZAF', 'SSD', 'SHN', 'SDN', 'SWZ', 'TZA', 'TGO', 'TUN', 'UGA', 'COD', 'ZMB', 'TZA', 'ZWE')
	# --------------------------------------------------------------------------------


	# -----------------------------------------------------
	# Load/prep data

	# load covariates
	edu = raster(eduFile)
	access = raster(accessFile)

	# load estimates
	vac = raster(vacFile)

	# load shapefile at country level
	map0 = shapefile(admin0File)

	# load shapefile at admin-1 level
	map1 = shapefile(admin1File)

	# load the ground cover data
	lakes = shapefile(shapeFileLakes)

	# subset to selected countries if specified above
	if (subset[1]!='All') { 
		map0 = map0[map0@data$'GID_0' %in% subset,]
		map1 = map1[map1@data$'GID_0' %in% subset,]
		edu = crop(edu, extent(map0))
		access = crop(access, extent(map0))
		vac = crop(vac, extent(map0))
		edu = mask(edu, map0)
		access = mask(access, map0)
		vac = mask(vac, map0)
	} 

	# crop shapefiles to vaccine coverage extent if all countries
	if (subset[1]=='All') { 
		map0 = crop(map0, extent(vac))
		map1 = crop(map1, extent(vac))
	}

	# project covariates to estimate pixels to ensure merge
	edu = projectRaster(edu, vac)
	access = projectRaster(access, vac)

	# crop to remove white space
	edu = crop(edu, extent(vac))
	access = crop(access, extent(vac))

	# mask out countries with no coverage estimates
	edu = mask(edu, vac)
	access = mask(access, vac)

	# mask out bodies of water
	edu = mask(edu, lakes, inverse=TRUE)
	access = mask(access, lakes, inverse=TRUE)
	vac = mask(vac, lakes, inverse=TRUE)

	# aggregate resolution if specified
	if (decrease_res) { 
		edu = aggregate(edu, fact=res_factor)
		access = aggregate(access, fact=res_factor)
		vac = aggregate(vac, fact=res_factor)
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
	# -----------------------------------------------------


	# ----------------------------------------------------------------------------------
	# Identify the admin1 for each pixel number

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
	# Simplify maps (after extraction)
	map0tmp = map0@data
	map0S = gSimplify(map0, tol=.01, topologyPreserve=TRUE)
	map0S = as(map0, 'SpatialPolygonsDataFrame')
	map0S@data = map0tmp
	map1tmp = map1@data
	map1S = gSimplify(map1, tol=.01, topologyPreserve=TRUE)
	map1S = as(map1S, 'SpatialPolygonsDataFrame')
	map1S@data = map1tmp
	# ------------------------------------------------------


	# -------------------------------
	# Run analyses
	# This script creates the following objects:
	# 'lmFit','lmeFit','explainedVariancesOverall', 'explainedVariancesAdmin1'
	source('./analyze_correlates.r')
	# -------------------------------


	# --------------------------------------------------------
	# Save
	objs = NULL
	for (o in c('access','edu','vac','data','map0S','map1S',
				'lmFit','lmeFit','explainedVariancesOverall',
				'explainedVariancesAdmin1')) {
		if (o %in% ls()) objs = c(objs, o)
	}
	save(list=objs, file=outFile)
	# -------------------------------------------------------

	
# ------------
# end function
}
# ------------
