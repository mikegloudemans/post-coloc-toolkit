require(rjson)
require(dplyr)

############################################################
### Filter colocalization results
############################################################

# Config parameters:
#
# "post_hoc_filter" :
# {
# 	"full_results_file": "results_file",	# Optional
# 	"out_results_file": "out_results_file",	# Optional
#
#	"kept_gwas": ["gwas1", "gwas2", ...],	# Optional ; CANNOT include both "kept_gwas" and "removed_gwas"
#	"kept_eqtl": ["eqtl1", "eqtl2", ...],	# Optional ; CANNOT include both "kept_eqtl" and "removed_eqtl"
#
#	"removed_gwas": ["gwas3"],		# Optional
#	"removed_eqtl": ["eqtl4"],		# Optional
#
#	"gwas_pval_threshold":			# Optional
#	{
#		"standard": 5e-8,		# Not optional if "gwas_pval_threshold" is included
#		"exceptions":			# Always optional
#		{
#			"special-gwas1": 1e-5,
#			"special-gwas2": 1e-5
#		}
#	},
#	"eqtl_pval_threshold":			# Optional
#	{
#		"standard": 1e-5,		# Not optional if "eqtl_pval_threshold" is included
#		"exceptions": {}		# Always optional
#	}
#
# }
# If included, "post_hoc_filter" is a JSON object that contains additional parameters
# that will be used for running this filtering step.
#
#	"full_results_file", "full_errors_file", and "full_skips_file" are file locations of
#	the results, errors, and skipped variants files, in the format produced by the "cat_results"
#	tool, relative to the top directory of this project.
#
#	"kept_gwas" is, optionally, a list of GWAS traits that should be included while all others are
#	filtered out. "removed_gwas", by contrast, is a list of GWAS traits that should be removed and
#	all others kept. To avoid conflicts, no more than one of "kept_gwas" and "removed_gwas" can be
#	specified in the options. "kept_eqtl" and "removed_eqtl" apply in a similar way. The format of the
#	GWAS and eQTL trait names in these fields should be identical to the strings used in the "eqtl_file"
#	and "trait" columns of the input files being filtered.
#
#	"gwas_pval_threshold" is, optionally, a p-value threshold at or above which results should be removed
#	from the results table. It must always include at least one field, "standard", which specifies a floating-
#	point number to use as the default threshold for filtering. Then, an "exceptions" field can optionally
# 	be included, to designate specific GWAS that are exceptions to the normal filtering threshold. In the example
#	above, "special-gwas1" and "special-gwas2" have been given exceptions to give them a less stringent GWAS cutoff
#	threshold (1e-5 instead of the standard 5e-8 that will be used for the rest of the loci). "eqtl_pval_threshold"
#	is exactly analogous to "gwas_pval_threshold".
#
#		A typical use case for this: perhaps you ran colocalization for 10 traits at all loci meeting a 
#		loose threshold of GWAS p-value < 1e-5. In your final presentation of the results, however, you'd
#		like to focus only on results meeting the typical genome-wide signficance threshold of GWAS
# 		p-value < 5e-8. But perhaps there are also 3 traits of particular interest, for which you're
#		willing to dig through a few extra, less-promising loci to ensure that you don't miss an interesting
#		ones. You set an exception for these 3 traits to include all loci passing the 1e-5 cutoff.
#
#		We note that there is not yet a rigorous statistical approach for controlling the false positive
# 		rate in these colocalization trials; therefore, care should be taken to follow up all identified
# 		colocalizations and to consider them within the context of other available multi-omic data. 
# 		In general, however, lowering the p-value cutoff detects more possible causal colocalizations, at
#		the cost of also turning up more likely false positives.
#
#	MAYBE also include other raw files here for filtering too while we're at it...



post_hoc_filter = function(config_file, input_file, output_file)
{
	config = fromJSON(file=config_file)$post_hoc_filter
	config$input_file = input_file
	config$output_file = output_file
	
	# Load results, errors, skips files
	results = load_post_hoc_filter_input_file(config)

	pre_results_dim = dim(results)[1]

	# If specified, filter results down to a limited set of GWAS and/or eQTL studies 

	if ("kept_gwas" %in% names(config))
	{
		results = filter_by_gwas(results, config$kept_gwas, keep=TRUE)
	}
	else if ("removed_gwas" %in%  names(config))
	{
		results = filter_by_gwas(results, config$removed_gwas, keep=FALSE)
	}

	if ("kept_eqtl" %in% names(config))
	{
		results = filter_by_eqtls(results, config$kept_eqtl, keep=TRUE)
	}
	else if ("removed_eqtl" %in%  names(config))
	{
		results = filter_by_gwas(results, config$removed_eqtl, keep=FALSE)
	}

	results = apply_pval_filter(results, config)

	# Display warning if not a single result was removed.
	if (dim(results)[1] == pre_results_dim)
	{
		print("Warning: No tests were removed during post-hoc filtering.")
	}

	# Write filtered output files
	write.table(results, config$output_file, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
}

# Remove rows from "data" that are not from one of the 
# designated GWAS.
filter_by_gwas = function(data, gwas, keep=TRUE)
{
	new_data = data[keep == (data$gwas_file %in% gwas),]
	return(new_data)
}

# Remove rows from "data" that are not from one of the 
# designated eQTL tissues.
filter_by_eqtl = function(data, eqtl, keep=TRUE)
{
	new_data = data[keep == (data$eqtl_file %in% eqtl),]
	return(new_data)
}

# create a function that tests GWAS pval / eQTL trait rows for validity
pval_passing = function(x, threshold_set)
{
	pvalue = as.numeric(x[1])
	trait = unlist(x[2])
	if (pvalue < threshold_set$standard)
	{
		return(TRUE)
	}
	else if (("exceptions" %in% names(threshold_set)) && (trait %in% names(threshold_set$exceptions)))
	{
		if (pvalue < threshold_set$exceptions[[trait]])
		{
			return(TRUE)
		}
	}
	return(FALSE)
}

load_post_hoc_filter_input_file = function(config)
{
	t = read.table(config$input_file, header=TRUE)

	t = t %>% select(ref_snp, eqtl_file, feature, n_snps, neg_log_gwas_pval, neg_log_eqtl_pval, gwas_file, score)

	return(t)	
}

# Function that tests GWAS pval / eQTL trait rows for validity
pval_passing = function(x, threshold_set)
{
	pvalue = 10^(-as.numeric(x[1]))
	trait = unlist(x[2])
	if (pvalue < threshold_set$standard)
	{
		return(TRUE)
	}
	else if (("exceptions" %in% names(threshold_set)) && (trait %in% names(threshold_set$exceptions)))
	{
		if (pvalue < threshold_set$exceptions[[trait]])
		{
			return(TRUE)
		}
	}
	return(FALSE)
}


apply_pval_filter = function(results, config)
{
	filtered_results = results
	if ("gwas_pval_threshold" %in% names(config))
	{
		filtered_results = filtered_results[apply(filtered_results[c("neg_log_gwas_pval", "gwas_file")], 1, FUN=pval_passing, threshold = config$gwas_pval_threshold),] 
	}
	if ("eqtl_pval_threshold" %in% names(config))
	{
		filtered_results = filtered_results[apply(filtered_results[c("neg_log_eqtl_pval", "eqtl_file")], 1, FUN=pval_passing, threshold = config$eqtl_pval_threshold),] 
	}

	return(filtered_results)
}

args = commandArgs(trailingOnly=TRUE)

config_file = args[1]
input_file = args[2]
output_file = args[3]
post_hoc_filter(config_file, input_file, output_file)
