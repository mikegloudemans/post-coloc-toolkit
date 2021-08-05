suppressWarnings(suppressMessages(require(rjson)))
suppressWarnings(suppressMessages(require(dplyr)))

############################################################
### Filter colocalization results
############################################################

classify_results = function(config_file, input_file, output_file, summary_file)
{
	# TODO: Have some sort of default way that this is done, if no custom file
	# is specified -- maybe just a default config file that is downloaded along with the other
	# ones

	# TODO: Make it so the config file for this can either be specified separately OR
	# defined directly in the config files

	# Load config file
	config = fromJSON(file=config_file)$classify_results
	config$input_file = input_file
	config$output_file = output_file
	config$summary_file = summary_file

	# Load results table
	results = load_results_file(config)
	
	# TODO: Validate rules too...

	# Apply rules, one at a time
	rule_list = config$rules
	for (rule_name in names(rule_list))
	{
		rule = rule_list[[rule_name]]
		# Add column tagging loci based on this rule...
		if (rule$type == "num_colocs")
		{
			results[[rule_name]] = class_by_num_coloc(results, rule, rule_name)
		} else if (rule$type == "specificity")
		{
			results[[rule_name]] = class_by_column_specificity(results, rule, rule_name)
		}

		# All loci should belong to a group at this point; if not, the groups are misspecified
		if(!(sum(results[[rule_name]] == "") == 0))
		{
			print(sprintf("Warning: not all loci in have been assigned to a group in %s.
				Check to make sure rules define the entire space of loci.", rule_name))
		}
		locus_classes = suppressMessages(suppressWarnings(results %>% group_by(!!as.name(rule_name)) %>% summarize(num_loci=length(unique(locus)))))
		suppressWarnings(write.table(locus_classes, file = gsub("_completion_indicator.tmp", sprintf("_%s.txt", rule_name), config$summary_file), sep="\t", quote=FALSE, row.names=FALSE,col.names=TRUE))
	}

	# Output SNP table with loci
	write.table(results, config$output_file, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

	system(sprintf("touch %s", config$summary_file))
}


# Perform a sort based on the number of candidate genes
# and the total number of colocalized genes
class_by_num_coloc = function(results, rule, rule_name)
{
	loci_list = unique(results$locus)

	class_membership = rep("", length(loci_list))

	# For each feature, test whether it passed colocalization threshold at this locus 
	summary = suppressMessages(suppressWarnings(results %>% group_by(locus, ensembl) %>% summarize(colocs=sum(coloc_status=="coloc"))))

	# Now summarize the number of colocalized genes and the number of candidate genes at each locus
	coloc_counts = suppressMessages(suppressWarnings(summary %>% group_by(locus) %>% summarize(num_coloc_genes = sum(colocs > 0), num_candidate_genes=length(ensembl))))

	# Each of these rule types takes as input a value X, a list of indices of candidate loci, and the
	# metadata about the number of colocalizations and tested genes at every locus.
	# It outputs the subset of these indices for which the corresponding loci pass the test.
	num_candidates_equals = function(x, loci, coloc_matrix)
	{
		return(loci[loci %in% (coloc_matrix[coloc_counts$num_candidate_genes == as.numeric(x),]$locus)])
	}
	num_colocs_equals = function(x, loci, coloc_matrix)
	{
		return(loci[loci %in% (coloc_matrix[coloc_counts$num_coloc_genes == as.numeric(x),]$locus)])
	}
	num_candidates_greater_than = function(x, loci, coloc_matrix)
	{
		return(loci[loci %in% (coloc_matrix[coloc_counts$num_candidate_genes > as.numeric(x),]$locus)])
	}
	num_colocs_greater_than = function(x, loci, coloc_matrix)
	{
		return(loci[loci %in% (coloc_matrix[coloc_counts$num_coloc_genes > as.numeric(x),]$locus)])
	}
	num_candidates_less_than = function(x, loci, coloc_matrix)
	{
		return(loci[loci %in% (coloc_matrix[coloc_counts$num_candidate_genes < as.numeric(x),]$locus)])
	}
	num_colocs_less_than = function(x, loci, coloc_matrix)
	{
		return(loci[loci %in% (coloc_matrix[coloc_counts$num_coloc_genes < as.numeric(x),]$locus)])
	}

	for (type in names(rule$categories))
	{
		# Keep track of all the loci passing each rule
		pass = 1:max(coloc_counts$locus)

		# Now go through all the rule types for this class to see which loci pass that rule
		if ("num_candidates_equals" %in% names(rule$categories[[type]]))
		{
			pass = num_candidates_equals(rule$categories[[type]]["num_candidates_equals"], pass, coloc_counts)
		}
		if ("num_colocs_equals" %in% names(rule$categories[[type]]))
		{
			pass = num_colocs_equals(rule$categories[[type]]["num_colocs_equals"], pass, coloc_counts)
		}
		if ("num_candidates_greater_than" %in% names(rule$categories[[type]]))
		{
			pass = num_candidates_greater_than(rule$categories[[type]]["num_candidates_greater_than"], pass, coloc_counts)
		}
		if ("num_colocs_greater_than" %in% names(rule$categories[[type]]))
		{
			pass = num_colocs_greater_than(rule$categories[[type]]["num_colocs_greater_than"], pass, coloc_counts)
		}
		if ("num_candidates_less_than" %in% names(rule$categories[[type]]))
		{
			pass = num_candidates_less_than(rule$categories[[type]]["num_candidates_less_than"], pass, coloc_counts)
		}
		if ("num_colocs_less_than" %in% names(rule$categories[[type]]))
		{
			pass = num_colocs_less_than(rule$categories[[type]]["num_colocs_less_than"], pass, coloc_counts)
		}


		# Make sure no locus has been double-classified; this would be a mistake
		if (!(sum(class_membership[which(loci_list %in% pass)] != "") == 0))
		{
			print(sprintf("Class specification error with rule %d. Some loci belong to more than one class.", rule_name))	
			# TODO: Give more details about which rules have collided.
			stop()
		}
		class_membership[which(loci_list %in% pass)] = type

	}

	all_classes = class_membership[match(results$locus, loci_list)]

	return(all_classes)
}

# Perform a sort based on the number of candidate genes
# and the total number of colocalized genes
class_by_column_specificity = function(results, rule, coloc_threshold, rule_name)
{
	loci_list = unique(results$locus)

	class_membership = rep("None", length(loci_list))

	# Figure out which tissues had strong, weak, no colocs at each locus
	tissue_coloc = suppressMessages(suppressWarnings(results %>% group_by(locus, !!as.name(rule$column)) %>% summarize(has_coloc = as.numeric(sum(coloc_status == "coloc") > 0))))

	all_colocs = tissue_coloc[tissue_coloc$has_coloc == TRUE,] 
	coloc_loci = unique(all_colocs$locus)

	# Run through the rules backwards, in ascending order of priority,
	# since some rules may satisfy more than one class
	for (class in rev(rule$categories))
	{
		# Keep track of all the loci passing each rule
		pass = sapply(loci_list, function(x)
	       	{
			this_locus = all_colocs %>% filter(locus == x)
			locus_tissues = unique(this_locus[[rule[["column"]]]])

			# If no colocalizations, it's just classed as None, the default
			if (length(locus_tissues) == 0)
			{
				return(FALSE)
			}

			if ("contains_exactly" %in% names(class))
			{
				if ((sum(!(class$contains_exactly %in% locus_tissues)) != 0) || (sum(!(locus_tissues %in% class$contains_exactly))!=0))
				{
					return(FALSE)
				}	
			}
			if ("contains_all" %in% names(class))
			{
				if (sum(!(class$contains_all %in% locus_tissues)) != 0)
				{
					return(FALSE)
				}	
			}
			if ("contains_none" %in% names(class))
			{
				if (sum(class$contains_none %in% locus_tissues) != 0)
				{
					return(FALSE)
				}	
			}
			if ("contains_some" %in% names(class))
			{
				if (sum(class$contains_some %in% locus_tissues) == 0)
				{
					return(FALSE)
				}	
			}
			if ("contains_only" %in% names(class))
			{
				if (sum(!(locus_tissues %in% class$contains_only)) != 0)
				{
					return(FALSE)
				}	
			}
			return(TRUE)
	        })
		
		# NOTE: Some previous designations may be overwritten, since we're applying the
		# rules in ascending priority order
		class_membership[pass] = class$class_name
	}		
	all_classes = class_membership[match(results$locus, loci_list)]
	return(all_classes)

}

load_results_file = function(config)
{
	d = read.table(file=config$input_file, header=TRUE, sep="\t")

	return(d)
}

args = commandArgs(trailingOnly=TRUE)

config_file = args[1]
input_file = args[2]
output_file = args[3]
summary_file = args[4]
classify_results(config_file, input_file, output_file, summary_file)
