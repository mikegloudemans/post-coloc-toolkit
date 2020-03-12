require(reshape2)
require(ggplot2)
require(dplyr)
require(readr)
require(rjson)

############################################################
### Concatenate colocalization results
############################################################

# Required cat_config parameters (if this step is not skipped):
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
# "skip_steps": ["cat_results", ...]
# If "cat_results" is not included in the "skip_steps" list
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

cat_results = function(config)
{
	# Create handle for config settings specific to this module`
	cat_config = config$tool_settings$cat_results 

	# User specifies a list of directories containing the required
	# results, R goes through the directories to find the relevant files.
	# We then load and concatenate them into one single table.
	results = get_cat_results_input(cat_config$input_dirs)
	results = munge_results(results)

	errors = get_error_logs(cat_config$input_dirs)
	skips = get_skip_logs(cat_config$input_dirs)

	# If specified, get plot-friendly names for each GWAS and eQTL tissue/trait.
	results$gwas_display_name = get_gwas_display_names(results$gwas_file, cat_config)
	results$eqtl_display_name = get_eqtl_display_names(results$eqtl_file, cat_config)

	# Output a single table each for results, errors, and skipped variants.
	write.table(results, file=paste(config$output_dir, cat_config$cat_results_coloc_out_file, sep="/"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	write.table(errors, file=paste(config$output_dir, cat_config$cat_results_errors_out_file, sep="/"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	write.table(skips, file=paste(config$output_dir, cat_config$cat_results_skips_out_file, sep="/"), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
}

get_cat_results_input = function(raw_coloc_output_dirs)
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
	for (d in raw_coloc_output_dirs)
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
	for (d in raw_coloc_output_dirs)
	{
		files = list.files(path=d, full.names = TRUE, recursive = TRUE)
		skipped_files = files[grepl("skipped_variants.txt", files)]
		
		for (sf in skipped_files)
		{
			skip_tab[[i]] = read.table(sf, header=FALSE, sep="\t", fill=TRUE, col.names = c("gwas_file", "eqtl_file", "snp.chrom", "snp.pos", "feature", "error", "gwas_data"))
			i = i+1
		}
	}
	skips = do.call(rbind, skip_tab)
	return(skips)
}

# Optionally, create "display names" for each GWAS that
# can be used later in plots and tables.
get_gwas_display_names = function(gwas_names, cat_config)
{
	gwas_display = as.character(gwas_names)
	if ("gwas_display_names" %in% names(cat_config))
	{
		for (full_name in names(cat_config$gwas_display_names))
		{
			gwas_display[gwas_display == full_name] = cat_config$gwas_display_names[[full_name]]
		}
	}
	return(gwas_display)
}

# Optionally, create "display names" for each eQTL that
# can be used later in plots and tables.
get_eqtl_display_names = function(eqtl_names, cat_config)
{
	eqtl_display = as.character(eqtl_names)
	if ("eqtl_display_names" %in% names(cat_config))
	{
		for (full_name in names(cat_config$eqtl_display_names))
		{
			eqtl_display[eqtl_display == full_name] = cat_config$eqtl_display_names[[full_name]]
		}
	}
	return(eqtl_display)
}

# A temporary function just for now, to make existing
# results headers compatible with the pipeline
munge_results = function(results)
{
	mr = results
	colnames(mr)[which(colnames(mr) == "base_gwas_file")] = "gwas_file"
	colnames(mr)[which(colnames(mr) == "X.log_gwas_pval")] = "gwas_neg_log_pvalue"
	colnames(mr)[which(colnames(mr) == "X.log_eqtl_pval")] = "eqtl_neg_log_pvalue"
	return(mr)
}
