require(rjson)
require(rjson)
require(dplyr)

summarize_filter_counts = function(config_file, input_file, output_file)
{
	config = fromJSON(file=config_file)$summarize_filter_counts
	config$input_file = input_file
	config$output_file = output_file

	# Get pre-filtering file, if needed
	# Annotate it with loci
	# Filter it down by post-hoc-filter criteria
	if ("pre_overlap_file" %in% names(config))
	{
		all_loci = 
	}

	# And get coloc file, which will show pre- and post-coloc

	# Output # of loci and # of genes passing at each level

	# Do it on an individual GWAS level too, again # loci and # genes

}

# Remove rows from "data" that are not from one of the 
# designated GWAS.
filter_by_gwas = function(data, gwas, keep=TRUE)
{
	new_data = data[keep == (data$gwas_trait %in% gwas),]
	return(new_data)
}

# Remove rows from "data" that are not from one of the 
# designated eQTL tissues.
filter_by_qtl = function(data, qtl, keep=TRUE)
{
	new_data = data[keep == (data$qtl_file %in% qtl),]
	return(new_data)
}

get_coloc_status = function(data, config)
{
	new_data = data
	new_data$coloc_status = "none"
	new_data$score = as.numeric(new_data$score)
	new_data$coloc_status[new_data$score > as.numeric(config$colocalization_threshold)] = "coloc"
	return(new_data)
}

load_post_hoc_filter_input_file = function(config)
{
	t = read.table(config$input_file, header=TRUE, stringsAsFactors=FALSE)

	t = t %>% select(ref_snp, qtl_file, feature, n_snps, neg_log_gwas_pval, neg_log_qtl_pval, gwas_trait, score, ensembl)

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
		filtered_results = filtered_results[apply(filtered_results[c("neg_log_gwas_pval", "gwas_trait")], 1, FUN=pval_passing, threshold = config$gwas_pval_threshold),] 
	}
	if ("qtl_pval_threshold" %in% names(config))
	{
		filtered_results = filtered_results[apply(filtered_results[c("neg_log_qtl_pval", "qtl_file")], 1, FUN=pval_passing, threshold = config$qtl_pval_threshold),] 
	}

	return(filtered_results)
}

apply_snp_count_filter = function(results, config)
{
	filtered_results = results
	if ("snp_count_min_threshold" %in% names(config))
	{
		filtered_results = filtered_results %>% filter(as.numeric(as.character(filtered_results$n_snps)) >= as.numeric(config$snp_count_min_threshold))
	}
	return(filtered_results)
}

args = commandArgs(trailingOnly=TRUE)

config_file = args[1]
input_file = args[2]
output_file = args[3]
summarize_filter_counts(config_file, input_file, output_file)


