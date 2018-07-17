# --------------------------------------------------------------------
# David Phillips
#
# 7/17/2018
# Loop over different runs and collate a table of explained variances
# --------------------------------------------------------------------


# ------------------
# Set up R
rm(list=ls())
library(data.table)
library(stringr)
# ------------------

# ---------------------------------------------------------------
# Directories and files
dir = '/home/j/temp/davidp6/vaccine_coverage_correlates/data/'

files = c('input_data_UGAZMBresolution1.rdata',
'input_data_UGAZMBresolution4.rdata',
'input_data_Allresolution16.rdata',
'input_data_Allresolution24.rdata',
'input_data_Africaresolution24.rdata',
'input_data_Africaresolution16.rdata',
'input_data_Allresolution1.rdata',
'input_data_UGAZMBresolution16.rdata',
'input_data_Africaresolution1.rdata',
'input_data_Allresolution4.rdata',
'input_data_Allresolution8.rdata',
'input_data_Africaresolution4.rdata',
'input_data_Africaresolution8.rdata',
'input_data_UGAZMBresolution8.rdata')

# output files
outFile1 = paste0(dir, 'evTableAdmin1.csv')
outFile2 = paste0(dir, 'evTableOverall.csv')
# ---------------------------------------------------------------


# --------------------------------------------------------------------
# Loop over files and append
i=1
for (f in files) { 
	print(paste('Loading file:', f))
	
	# load data
	load(paste0(dir, f))
	
	# parse name
	cleaned = gsub('input_data_', '', f)
	cleaned = gsub('.rdata', '', cleaned)
	parsed = str_split(cleaned, 'resolution',)[[1]]
	geo = parsed[1]
	res = parsed[2]

	# label runs
	explainedVariancesAdmin1[, geography:=geo]
	explainedVariancesAdmin1[, resolution:=res]
	explainedVariancesOverall[, geography:=geo]
	explainedVariancesOverall[, resolution:=res]
	
	if (i==1) evTableA = explainedVariancesAdmin1
	if (i==1) evTableO = explainedVariancesOverall
	if (i>1) evTableA = rbind(evTableA, explainedVariancesAdmin1)
	if (i>1) evTableO = rbind(evTableO, explainedVariancesOverall)
	
	i=i+1
}

# reshape overall table
evTableO = dcast.data.table(evTableO, geography+resolution~variable, value.var='explained_variance')
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# Print to screen
means = tmp[, lapply(.SD, mean), by=c('iso3','geography','resolution'),
			.SDcols=c('edu','access','supply')]

uga = means[iso3=='UGA']
uga[, dist:=sqrt((edu-.36)^2 + (access-.2)^2 + (supply-(.25+.18))^2)]
uga[, dist1:=sqrt((edu-.36)^2)]
uga[, dist2:=sqrt((access-.2)^2)]
uga[, dist3:=sqrt((supply-(.25+.18))^2)]

zmb = means[iso3=='ZMB']
zmb[, dist:=sqrt((edu-.4)^2 + (access-.22)^2 + (supply-(.15+.23))^2)]
zmb[, dist1:=sqrt((edu-.4)^2)]
zmb[, dist2:=sqrt((access-.22)^2)]
zmb[, dist3:=sqrt((supply-(.15+.23))^2)]
# --------------------------------------------------------------------


# ---------------------------------------------
# Save
write.csv(evTableA, outFile1, row.names=FALSE)
write.csv(evTableO, outFile2, row.names=FALSE)
# ---------------------------------------------
