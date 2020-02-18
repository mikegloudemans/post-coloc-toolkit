require(reshape2)
require(ggplot2)
require(dplyr)
require(readr)
require(rjson)

############################################################
### Concatenate colocalization results
############################################################

# Required config parameters (if this step is not skipped):
#
#
#
# "raw_coloc_output_dirs": ["dir1", "dir2", "dir3", ...]
# Paths to directories containing all raw coloc output files.
# Can be absolute or relative paths; relative paths must be
# specified relative to the top-level directory of this project
# containing the master workflow script `run_full_workflow.py`.
# Specified directories need not contain only colocalization
# results files; can contain other files too. Specified folders will be
# explored recursively and all children files containing colocalization
# results will be included.
#
#
#
# "output_dir": "dir4"
# Path to the directory where all output files should be placed, including
# the concatenated colocalization results produced by this script
#
#
#


# Optional config parameters:
#
#
#
# "skip_steps": ["cat-results", ...]
# If "cat-results" is not included in the "skip_steps" list
# parameter, then it will run by default. If included,
# then this step will be skipped.
#
#
#
# "gwas_display_names":
# {
# 	"gwas1_long_name": "gwas1_pretty_name",
#	"gwas2_long_name": "gwas2_pretty_name",
#	...
# } 
# If specified, "gwas_display_names" is a JSON object mapping
# longer GWAS file names to names suitable for plotting and
# publication, which may include spaces or any characters that
# are part of a valid R string. Not all GWAS need to have a mapping;
# a subset or even none can be included. If none are included, the
# gwas_display_name column in the output table will simply be the
# same as the gwas_file column.
#
#
#
# "eqtl_display_names":
# {
# 	"eqtl1_long_name": "eqtl1_pretty_name",
#	"eqtl2_long_name": "eqtl2_pretty_name",
# 	...
# } 
# If specified, "eqtl_display_names" is a JSON object mapping
# longer eQTL file names to names suitable for plotting and
# publication, which may include spaces or any characters that
# are part of a valid R string. Not all eQTL need to have a mapping;
# a subset or even none can be included. If none are included, the
# eqtl_display_name column in the output table will simply be the
# same as the eqtl_file column.
#
#
#


main = function()
{
	# Load config file
	config_file = commandArgs(trailingOnly=TRUE)[1]
	config = validate_config(config_file)

	# User specifies a list of directories containing the required
	# results, R goes through the directories to find the relevant files.
	# We then load and concatenate them into one single table.
	results = get_raw_coloc_results(config$raw_coloc_output_dirs)
	errors = get_error_logs(config$raw_coloc_output_dirs)
	skips = get_skip_logs(config$raw_coloc_output_dirs)

	# If specified, get plot-friendly names for each GWAS and eQTL tissue/trait.
	results$gwas_display_name = get_gwas_display_names(results$gwas_file, config)
	results$eqtl_display_name = get_eqtl_display_names(results$eqtl_file, config)

	# Output a single table each for results, errors, and skipped variants.
	write.table(results, file=paste(config$output_dir, "raw_coloc_table.txt", sep="/"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	write.table(errors, file=paste(config$output_dir, "raw_coloc_errors_table.txt", sep="/"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	write.table(skips, file=paste(config$output_dir, "raw_coloc_skips_table.txt", sep="/"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
}

validate_config = function(config_file)
{
	# Load config file, specified as a command line parameter
	config = fromJSON(file=config_file)

	# Validate config file to be sure required parameters are present
	if (!("raw_coloc_output_dirs" %in% names(config)))
	{
		stop("Config ERROR: You must specify 'raw_coloc_output_dirs' in config file, 
		     or else skip the cat-results step.")
	}
	if (!("output_dir" %in% names(config)))
	{
		stop("Config ERROR: You must specify 'output_dir' in config file")
	}

	return(config)
}

get_raw_coloc_results = function(raw_coloc_output_dirs)
{
	tabs = list()
	i = 1
	for (d in raw_coloc_output_dirs)
	{
		files = list.files(path=d, full.names = TRUE, recursive = TRUE)

		# TODO: Expand this to allow it to include more results types than just FINEMAP files
		results_files = files[grepl("finemap_clpp_status.txt", files)]

		for (rf in results_files)
		{
			tabs[[i]] = read.table(rf, header=TRUE)
			i = i+1
		}	
	}
	# Combine results to a single table
	results = do.call(rbind, tabs)
	# It's possible we'll have exact duplicates if the traits weren't run all at once.
	# Remove them.
	results = results[!duplicated(results),]
	return(results)
}

get_error_logs = function(raw_coloc_output_dirs)
{
	err_tab = list()
	i = 1
	for (d in config$raw_coloc_output_dirs)
	{
		files = list.files(path=d, full.names = TRUE, recursive = TRUE)
		error_files = files[grepl("ERROR_variants.txt", files)]
		
		for (ef in error_files)
		{
			err_tab[[i]] = read.table(ef, header=FALSE, sep="\t", fill=TRUE, col.names = c("gwas_file", "eqtl_file", "snp.chrom", "snp.pos", "restrict_gene", "trait", "error", "error2"))
			err_tab[[i]]$error = paste(err_tab[[i]]$error, err_tab[[i]]$error2, sep="|")
			i = i+1
		}
	}
	errors = do.call(rbind, err_tab)
	return(errors)
}

get_skip_logs = function(raw_coloc_output_dirs)
{
	skip_tab = list()
	i = 1
	for (d in config$raw_coloc_output_dirs)
	{
		files = list.files(path=d, full.names = TRUE, recursive = TRUE)
		skip_files = files[grepl("skipped_variants.txt", files)]
		
		for (sf in skipped_files)
		{
			skip_tab[[i]] = read.table(paste(folder, skip, sep="/"), header=FALSE, sep="\t", fill=TRUE, col.names = c("gwas_file", "eqtl_file", "snp.chrom", "snp.pos", "feature", "error", "gwas_data"))
			i = i+1
		}
	}
	skips = do.call(rbind, skip_tab)
	return(skips)
}

# Optionally, create "display names" for each GWAS that
# can be used later in plots and tables.
get_gwas_display_names(gwas_names, config)
{
	gwas_names = as.character(gwas_names)
	gwas_display = gwas_names
	if ("gwas_display_names" %in% names(config))
	{
		for (full_name in names(config$gwas_display_names))
		{
			gwas_display[gwas_names == full_name] = config$gwas_display_names[[full_name]]
		}
	}
	return(gwas_display)
}

# Optionally, create "display names" for each eQTL that
# can be used later in plots and tables.
get_eqtl_display_names(eqtl_names, config)
{
	eqtl_names = as.character(eqtl_names)
	eqtl_display = eqtl_names
	if ("eqtl_display_names" %in% names(config))
	{
		for (full_name in names(config$eqtl_display_names))
		{
			eqtl_display[eqtl_names == full_name] = config$eqtl_display_names[[full_name]]
		}
	}
}

main()
