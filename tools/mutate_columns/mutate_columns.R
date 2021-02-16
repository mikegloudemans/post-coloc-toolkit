require(rjson)
require(dplyr)

############################################################
### Create / mutate new column names for display
############################################################

# Config parameters:
# TODO: 

mutate_columns = function(config_file, input_file, output_file)
{
	config = fromJSON(file=config_file)$mutate_columns
	config$input_file = input_file
	config$output_file = output_file
	
	# Load results, errors, skips files
	results = load_results_file(config)

	for (mutation in config$mutations)
	{
		column = results[[mutation[["in"]]]]
		results[[mutation[["out"]]]] = sapply(column, function(x)
			 {
				 if (x %in% names(mutation$map))
				 {
					 return(mutation$map[[x]])
				 }
				 else
				 {
					 return("NA")
				 }
			 })
	}

	# Write filtered output files
	write.table(results, config$output_file, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
}

load_results_file = function(config)
{
	t = read.table(config$input_file, header=TRUE, sep="\t")

	return(t)	
}

args = commandArgs(trailingOnly=TRUE)

config_file = args[1]
input_file = args[2]
output_file = args[3]
mutate_columns(config_file, input_file, output_file)
