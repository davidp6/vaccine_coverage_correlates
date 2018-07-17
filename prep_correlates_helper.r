# --------------------------------------------------------------------------------------------
# David Phillips
#
# 7/16/2018
# Helper script to run prep_correlates.r with command-line arguments (for running in parallel)
# Expected command-line arguments (see prep_correlates.r):
# Aguments 1-3: system variables
# Argument 4: decrease_res
# Argument 5: res_factor
# Argument 6-9: subset
# The current working directory should be the root of this repo
# --------------------------------------------------------------------------------------------


# ------------------------------------------
# Handle command line arguments
decrease_res = as.logical(commandArgs()[4])
res_factor = as.numeric(commandArgs()[5])
subset = commandArgs()[6]

# in case we pass a vector as the last arg
n = length(commandArgs())
if (n>6) { 
	for(i in seq(from=7, to=n)) {
		subset = c(subset, commandArgs()[i])
	}
}
# ------------------------------------------


# --------------------------------------------------------
# Run it

# print to screen
print(paste('Running:', decrease_res, res_factor, subset))

# load the prepCorrelates function
source('./prep_correlates.r')

# run the function with the command arguments
prepCorrelates(decrease_res, res_factor, subset)
# --------------------------------------------------------
