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
outFile1 = paste0(dir, 'evTableAdmin1.csv'))
outFile2 = paste0(dir, 'evTableOverall.csv')
# ---------------------------------------------------------------


# --------------------------------------------------------------------
# Loop over files and append
i=1
for (f in files) { 
	print(paste('Loading file:', f))
	
	# load data
	load(paste0(dir, f))
	
	if (i==1) evTableA = explainedVariancesAdmin1
	if (i>1) evTableA = rbind(evTable, explainedVariancesAdmin1)
	
	if (i==1) evTableO = explainedVariancesOverall
	if (i>1) evTableO = rbind(evTable, explainedVariancesOverall)

	# parse name
	cleaned = gsub('input_data_', '', f)
	cleaned = gsub('.rdata', '', cleaned)
	parsed = str_split(cleaned, 'resolution',)[[1]]
	geo = parsed[1]
	res = parsed[2]

	# label runs
	evTableA[, geography:=geo]
	evTableA[, resolution:=res]
	evTableO[, geography:=geo]
	evTableO[, resolution:=res]
	
	i=i+1
}
# --------------------------------------------------------------------


# ---------------------------------------------
# Save
write.csv(evTableA, outFile1, row.names=FALSE)
write.csv(evTableO, outFile2, row.names=FALSE)
# ---------------------------------------------
