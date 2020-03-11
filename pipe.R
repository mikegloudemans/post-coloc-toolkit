################################################
################################################
### General pipeline for post-coloc-toolkit
################################################
################################################

# Load all the tools
source("tools/cat_results/cat_results.R")
source("tools/validate_config/validate_config.R")

# Load config file as JSON object, from the specified input file

if(length(commandArgs(trailing=TRUE)) < 1)
{
	# Display an informative error if no config
	# file argument was specified
	stop("You must specify a JSON config file as a command-line argument.")
}

config_file = commandArgs(trailing=TRUE)[1]
config = fromJSON(file=config_file)

# Make sure config file is valid; this will throw an
# error and stop the program if it isn't
validate_config(config)

# Make the output directory if it doesn't already exist
if (!dir.exists(config$output_dir))
{
	dir.create(config$output_dir, recursive=TRUE)
	print(paste0("Created output directory ", config$output_dir))
}

#########################################
# Dispatch all the other scripts
# as requested
#########################################

# Concatenate results, if requested
if (!("skip_steps" %in% names(config)) || !("cat_results" %in% names(config$skip_steps))) 
{
	cat_results(config)
}


