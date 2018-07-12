# -----------------------------------------------------------------------------------------------
# David Phillips
#
# 7/3/2018
# Explore covariates and their relationship with vaccine coverage at the pixel level
# Note: cowplot may not be available for unix machines (it's used here to create complex legends)
# -----------------------------------------------------------------------------------------------


# --------------------
# Set up R
rm(list=ls())
library(data.table)
library(raster)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
# --------------------


# --------------------------------------------------------------
# Files and directories

j = ifelse(Sys.info()[1]=='Windows', 'J:', '/home/j')

# directory for output
dir = paste0(j, '/temp/davidp6/vaccine_coverage_correlates/')

# input file
inFile = paste0(dir, 'input_data_Africa.rdata')

# second input file if explained variances come from a different model run
# this loads first and assumes there are no model estimates in inFile1
inFile2 = paste0(dir, 'input_data_Africa_factor8.rdata')

# output file
graphFile = paste0(dir, 'vaccine_coverage_correlates_Africa.pdf')

# countries to subset to (ISO3 code(s), or "All")
subset = c('UGA','ZMB')
# --------------------------------------------------------------


# ------------------------------------------------
# Load all data
load(file=inFile2)
load(file=inFile)

# subset if specified
if (subset[1]!='All') { 
	data = data[iso3 %in% subset]
	map0S = map0S[map0S@data$'GID_0' %in% subset,]
	map1S = map1S[map1S@data$'GID_0' %in% subset,]
}
# ------------------------------------------------


# ---------------------------------------------------------------------------
# Clean data for graphing

# omit na to reduce extent
data = na.omit(data)

# limit access to 5 minutes and 10 days
data[, access_trimmed:=access2_mean_synoptic]
data[access2_mean_synoptic>14400, access_trimmed:=14400]
data[access2_mean_synoptic<5, access_trimmed:=5]

# bin pixels
data[, edu_bin:=cut(edu_mean_stage_2_mean_1y_2016_00_00, 
	quantile(edu_mean_stage_2_mean_1y_2016_00_00, c(0,.33,.66,1), na.rm=TRUE), 
	include.lowest=TRUE, labels=paste(c('Low','Medium','High'), 'Education'))]
data[, vac_bin:=cut(dpt3_cov_mean_raked_raster, 
	quantile(dpt3_cov_mean_raked_raster, c(0,.33,.66,1), na.rm=TRUE), 
	include.lowest=TRUE, labels=paste(c('Low','Medium','High'), 'Coverage'))]
data[, access_bin:=cut(access_trimmed, 
	quantile(access_trimmed, c(0,.33,.66,1), na.rm=TRUE), 
	include.lowest=TRUE, labels=paste(c('Low','Medium','High'), 'Access'))]

# make bivariate bins
data[, edu_vac_bin:=paste0(vac_bin, '/', edu_bin)]
data[, access_vac_bin:=paste0(vac_bin, '/', access_bin)]

# simplify maps
# map0tmp = map0@data
# map0 = gSimplify(map0, tol=.1, topologyPreserve=TRUE)
# map0 = as(map0, 'SpatialPolygonsDataFrame')
# map0@data = map0tmp
# map1tmp = map1@data
# map1S = gSimplify(map1, tol=.1, topologyPreserve=TRUE)
# map1S = as(map1S, 'SpatialPolygonsDataFrame')
# map1S@data = map1tmp

# fortify the map for ggplot
map0F = data.table(fortify(map0S, region='GID_0'))
map0F[, iso3:=id]
map1F = data.table(fortify(map1S))
map1F = merge(map1F, data.table(id=unique(map1F$id), admin1_id=map1S@data$GID_1), by='id')
map1F[, iso3:=str_sub(admin1_id,1,3)]

# # crop the map down to just the extent with data
# mapC = crop(map0, extent(vac))
# mapF = fortify(mapC)

# merge random effects to fortified map1
map1F = merge(map1F, explainedVariancesAdmin1, by='admin1_id')
# ---------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------
# Set up to graph

# define color order
colOrder1 = c('Low Coverage/Low Education', 'Low Coverage/Medium Education', 'Low Coverage/High Education',
			'Medium Coverage/Low Education', 'Medium Coverage/Medium Education', 'Medium Coverage/High Education',
			'High Coverage/Low Education', 'High Coverage/Medium Education', 'High Coverage/High Education')
colOrder2 = c('Low Coverage/Low Access', 'Low Coverage/Medium Access', 'Low Coverage/High Access',
			'Medium Coverage/Low Access', 'Medium Coverage/Medium Access', 'Medium Coverage/High Access',
			'High Coverage/Low Access', 'High Coverage/Medium Access', 'High Coverage/High Access')

# bivariate color scales
# GnBl from http://www.joshuastevens.net/cartography/make-a-bivariate-choropleth-map/
# BlRd, BrBl from https://www.slideshare.net/aileenbuckley/arc-gis-bivariate-mapping-tools-28903069
colsBlGn = c('#e8e8e8','#b8d6be','#73ae80','#b5c0da','#90b2b3','#5a9178','#6c83b5','#567994','#2a5a5b')
colstwRdBl = c('#FFFFC1','#E88E6F','#D5191B','#AAD8E8','#A1777F','#AD2E3F','#2979B6','#57628B','#804968')			
colstwBrBl = c('#E3EAD7','#EADB8F','#F3B201','#B4D1DE','#AEB3AD','#B06402','#4A9EC0','#376488','#344E4D')
colstwGnRd = c('#FFFFC1','#8DC983','#1E9041','#E98C6E','#988251','#487738','#D61A1B','#A63724','#75552C')

# select one of the above scales for education vs coverage
colGrid1 = colstwBrBl
names(colGrid1) = colOrder1

# select one of the above scales for access vs coverage
colGrid2 = colstwRdBl
names(colGrid2) = colOrder2

# extract univariate color scales from bivariate
cols1 = colGrid2[c(1,4,7)]
cols2 = colGrid1[c(1,2,3)]
cols3 = colGrid2[c(1,2,3)]

# overwrite univariate color scales
cols1 = brewer.pal(6, 'YlGnBu')
cols2 = brewer.pal(6, 'YlOrBr')
cols3 = brewer.pal(6, 'Reds')
cols4 = brewer.pal(6, 'Reds')
cols5 = brewer.pal(6, 'Reds')

# other colors
border = 'grey35'
# ----------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------
# Store graphs

# coverage estimates
p1 = ggplot(data, aes(y=y,x=x,fill=dpt3_cov_mean_raked_raster)) + 
	geom_tile() + 
	geom_path(data=map0F, aes(x=long, y=lat, group=group)
		, color=border, size=.025, inherit.aes=FALSE) + 
	scale_fill_gradientn('DPT3 Coverage', colors=cols1, na.value='white') + 
	labs(title='2016 DPT Coverage (Third Dose)') + 
	theme_void() + 
	theme(plot.title=element_text(hjust=1, face='plain'), 
		plot.subtitle=element_text(hjust=1), plot.margin=unit(c(2,0,2,0), 'cm'))
if (subset[1]!='All') p1 = p1 + facet_wrap(~iso3, scales='free')

		
# education
p2 = ggplot(data, aes(y=y,x=x,fill=edu_mean_stage_2_mean_1y_2016_00_00)) + 
	geom_tile() + 
	geom_path(data=map0F, aes(x=long, y=lat, group=group)
		, color=border, size=.05, inherit.aes=FALSE) + 
	scale_fill_gradientn('Mean Years\nof Education', colors=cols2, na.value='white') + 
	labs(title='2016 Maternal Educational Attainment') + 
	theme_void() + 
	theme(plot.title=element_text(hjust=1, face='plain'), 
		plot.subtitle=element_text(hjust=1), plot.margin=unit(c(2,0,2,0), 'cm'))
if (subset[1]!='All') p2 = p2 + facet_wrap(~iso3, scales='free')
	
# access	
p3 = ggplot(data, aes(y=y,x=x,fill=access_trimmed)) + 
	geom_tile() + 
	geom_path(data=map0F, aes(x=long, y=lat, group=group)
		, color=border, size=.025, inherit.aes=FALSE) + 
	scale_fill_gradientn('Travel Time', colors=cols3, trans='log', na.value='white', 
		breaks=c(10,60,360,1440,4320), labels=c('10 min', '1 hour', '6 hour', '1 day', '3 days')) + 
	labs(title='2015 Access to Settlements >50,000 Population') + 
	theme_void() + 
	theme(plot.title=element_text(hjust=1, face='plain'), 
		plot.subtitle=element_text(hjust=1), plot.margin=unit(c(2,0,2,0), 'cm'))
if (subset[1]!='All') p3 = p3 + facet_wrap(~iso3, scales='free')

# bivariate map of education vs coverage
p4a= ggplot(data, aes(y=y,x=x,fill=edu_vac_bin)) + 
	geom_tile() + 
	geom_path(data=map0F, aes(x=long, y=lat, group=group)
		, color=border, size=.025, inherit.aes=FALSE) + 
	scale_fill_manual('',values=colGrid1) + 
	labs(title='2016 DPT3 Coverage and Maternal Education') + 
	theme_void() + 
	theme(legend.position='none', plot.title=element_text(hjust=0, face='plain'), 
		plot.subtitle=element_text(hjust=0), plot.margin=unit(c(2,0,2,0), 'cm'))
if (subset[1]!='All') p4a = p4a + facet_wrap(~iso3, scales='free')

legendData = data.table(labs = names(colGrid1))
legendData[grepl('Low Education',labs), edu:=1]
legendData[grepl('Medium Education',labs), edu:=2]
legendData[grepl('High Education',labs), edu:=3]
legendData[grepl('Low Coverage',labs), vac:=1]
legendData[grepl('Medium Coverage',labs), vac:=2]
legendData[grepl('High Coverage',labs), vac:=3]
p4legend = ggplot(legendData, aes(y=vac, x=edu, fill=labs)) + 
	geom_tile() + 
	scale_fill_manual('', values=colGrid1) + 
	scale_x_continuous(breaks=c(1,2,3),labels=c('Low','Med.','High')) + 
	scale_y_continuous(breaks=c(1,2,3),labels=c('Low','Med.','High')) + 
	labs(y='Coverage',x='Education') + 
	theme_minimal(base_size=13) + 
	theme(legend.position='none', panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank())

p4 = ggdraw() + draw_plot(p4a) + draw_plot(p4legend, x=.5, y=.55, width=.2, height=.24)

# bivariate map of access vs coverage
p5a= ggplot(data, aes(y=y,x=x,fill=access_vac_bin)) + 
	geom_tile() + 
	geom_path(data=map0F, aes(x=long, y=lat, group=group)
		, color=border, size=.025, inherit.aes=FALSE) + 
	scale_fill_manual('',values=colGrid2) + 
	labs(title='2016 DPT3 Coverage and Access to Settlements') + 
	theme_void() + 
	theme(legend.position='none', plot.title=element_text(hjust=0, face='plain'), 
		plot.subtitle=element_text(hjust=0), plot.margin=unit(c(2,0,2,0), 'cm'))
if (subset[1]!='All') p5a = p5a + facet_wrap(~iso3, scales='free')

legendData = data.table(labs = names(colGrid2))
legendData[grepl('Low Access',labs), access:=1]
legendData[grepl('Medium Access',labs), access:=2]
legendData[grepl('High Access',labs), access:=3]
legendData[grepl('Low Coverage',labs), vac:=1]
legendData[grepl('Medium Coverage',labs), vac:=2]
legendData[grepl('High Coverage',labs), vac:=3]
p5legend = ggplot(legendData, aes(y=vac, x=access, fill=labs)) + 
	geom_tile() + 
	scale_fill_manual('', values=colGrid2) + 
	scale_x_continuous(breaks=c(1,2,3),labels=c('High','Med.','Low')) + 
	scale_y_continuous(breaks=c(1,2,3),labels=c('Low','Med.','High')) + 
	labs(y='Coverage',x='Access') + 
	theme_minimal(base_size=13) + 
	theme(legend.position='none', panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank())

p5 = ggdraw() + draw_plot(p5a) + draw_plot(p5legend, x=.5, y=.55, width=.2, height=.24)

# admin1-level correlation maps
p6 = ggplot(map1F, aes(x=long, y=lat, group=group, fill=edu)) + 
	geom_polygon() + 
	geom_path(color=border, size=.01) + 
	scale_fill_gradientn('% Explained by Demand', colors=cols2, na.value='white') + 
	labs(title='Percentage of 2016 DPT3 Coverage Explained by Demand', 
		subtitle='At Second Administrative Level') + 
	coord_fixed(ratio=1) + 
	theme_void() + 
	theme(plot.title=element_text(hjust=0, face='plain'), 
		plot.subtitle=element_text(hjust=0))
if (subset[1]!='All') p6 = p6 + facet_wrap(~iso3, scales='free')
	
p7 = ggplot(map1F, aes(x=long, y=lat, group=group, fill=supply)) + 
	geom_polygon() + 
	geom_path(color=border, size=.01) + 
	scale_fill_gradientn('% Explained by Supply', colors=cols1, na.value='white') + 
	labs(title='Percentage of 2016 DPT3 Coverage Explained by Supply', 
		subtitle='At Second Administrative Level') + 
	coord_fixed(ratio=1) + 
	theme_void() + 
	theme(plot.title=element_text(hjust=0, face='plain'), 
		plot.subtitle=element_text(hjust=0))
if (subset[1]!='All') p7 = p7 + facet_wrap(~iso3, scales='free')
# --------------------------------------------------------------------------------------------


# --------------------------------
# Save graphs
pdf(graphFile, height=5.5, width=10)
p1
p2
p3
p4
p5
p6
p7
dev.off()
# --------------------------------
