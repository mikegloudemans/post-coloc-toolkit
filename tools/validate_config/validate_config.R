# TODO: Specify all constants at the top, here, instead of making the
# coder hunt through the code for them

# TODO: Break validate_config file down into separate functions

validate_config = function(config)
{ 
	new_config = config

	# Verify that all essential parameters are present and valid
	if (!("output_dir" %in% names(config)))
	{
		stop("Config ERROR: You must specify 'output_dir' in config file")
	}
	if (!("tool_settings" %in% names(config)))
	{
		stop("Config ERROR: Config file 'tool_settings' must contain a 'tool_settings' object")
	}

	# Check all the individual components first

	# Note: skip_steps might not even be defined, but that's OK.

	# Check for cat_results
	if(!("cat_results" %in% names(config$skip_steps)))
	{
		if (!("cat_results") %in% names(config$tool_settings))
		{
			stop("Config ERROR: You must specify 'cat_results' in the 'tool_settings' object, 
			     or else skip the cat_results step.")
		}
			
		cat_config = config$tool_settings$cat_results
		if (!("input_dirs" %in% names(cat_config)))
		{
			stop("Config ERROR: You must specify 'cat_results_input_dirs' in the 'cat_results' object, 
			     or else skip the cat_results step.")
		}

		if (!("cat_results_coloc_out_file" %in% names(cat_config)))
		{
			cat_config$cat_results_coloc_out_file = "cat_results_coloc_out.txt"
		}
		if (!("cat_results_errors_out_file" %in% names(cat_config)))
		{
			cat_config$cat_results_errors_out_file = "cat_results_errors_out.txt"
		}
		if (!("cat_results_skips_out_file" %in% names(cat_config)))
		{
			cat_config$cat_results_skips_out_file = "cat_results_skips_out.txt"
		}
		new_config$tool_settings$cat_results = cat_config
	}
	
	# Check for post_hoc_filter
	if(!("post_hoc_filter" %in% names(config$skip_steps)))
	{
		if (!("post_hoc_filter") %in% names(config$tool_settings))
		{
			stop("Config ERROR: You must specify 'post_hoc_filter' in the 'tool_settings' object, 
			     or else skip the post_hoc_filter step.")
		}
		filter_config = config$tool_settings$post_hoc_filter

		if (("kept_gwas" %in% names(filter_config)) && ("removed_gwas" %in%  names(filter_config)))
		{
			stop("config error: you can't specify both 'kept_gwas' and 'removed_gwas' in the 'post_hoc_filter'
			     tool settings.")
		}	
		if (("kept_eqtl" %in% names(filter_config)) && ("removed_eqtl" %in%  names(filter_config)))
		{
			stop("config error: you can't specify both 'kept_eqtl' and 'removed_eqtl' in the 'post_hoc_filter'
			     tool settings.")
		}

		if ((("gwas_pval_threshold" %in% names(filter_config)) && !("standard" %in% names(filter_config$gwas_pval_threshold))))
		{
			stop("config error: must include a 'standard' gwas threshold if 'gwas_pval_threshold' is specified
			     in the 'post_hoc_filter' tool settings.")
		}
		if ((("eqtl_pval_threshold" %in% names(filter_config)) && !("standard" %in% names(filter_config$eqtl_pval_threshold))))
		{
			stop("config error: must include a 'standard' eqtl threshold if 'eqtl_pval_threshold' is specified
			     in the 'post_hoc_filter' tool settings.")
		}

		# If no output file specified, use a default one
		if (!("out_file" %in% names(filter_config)))
		{
			filter_config$out_file = "post_hoc_filter_out.txt"
		}
		new_config$tool_settings$post_hoc_filter = filter_config
	}
	
	# Check for add_hgnc_names
	if(!("add_hgnc_names" %in% names(config$skip_steps)))
	{
		if (!("add_hgnc_names") %in% names(config$tool_settings))
		{
			stop("Config ERROR: You must specify 'add_hgnc_names' in the 'tool_settings' object, 
			     or else skip the add_hgnc_names step.")
		}

		hgnc_config = config$tool_settings$add_hgnc_names

		# Figure out the input ensembl to hgnc file
		if (("ensembl_to_hgnc_map_file" %in% names(hgnc_config)))
		{
			if (!(("ensembl_col_index" %in% names(hgnc_config)) && ("hgnc_col_index" %in% names(hgnc_config))))
			{
				stop("Config ERROR: When using a custom Ensembl-to-HGNC mapping, you
				     must specify 'ensembl_col_index' and 'hgnc_col_index' in the 'add_hgnc_names' parameters.")
			}
		}
		else
		{
			# Use default file if no specific mapping provided
			hgnc_config$ensemble_to_hgnc_map_file = "data/hgnc/ensembl_to_hgnc.txt"
			hgnc_config$ensembl_col_index = 1
			hgnc_config$hgnc_col_index = 3
		}
		
		# If no output file specified, use a default one
		if (!("out_file" %in% names(hgnc_config)))
		{
			hgnc_config$out_file = "add_hgnc_names_out.txt"
		}

		new_config$tool_settings$add_hgnc_names = hgnc_config
	}

	# Check for add_rsids
	if(!("add_rsids" %in% names(config$skip_steps)))
        {
		if (!("add_rsids") %in% names(config$tool_settings))
		{
			stop("Config ERROR: You must specify 'add_rsids' in the 'tool_settings' object, 
			     or else skip the add_rsids step.")
		}

		rsid_config = config$tool_settings$add_rsids
		
		# Make sure the columns are marked if a custom rsid file is being used.
		if (("dbsnp_file" %in% names(rsid_config)))
		{
			if (!("rsid_col_index" %in% names(rsid_config)))
			{
				stop("Config ERROR: When using a custom rsid mapping for 'add_rsids', you
				     must specify 'rsid_col_index' in the config file.")
			}
			if (!(("rsid_col_index" %in% names(rsid_config)) && ("chr_col_index" %in% names(rsid_config)
										  && ("pos_col_index" %in% names(rsid_config)))))
			{
				if (("use_tabix" %in% names(rsid_config)) && (rsid_config$use_tabix == "FALSE"))
				stop("Config ERROR: When using a custom rsid mapping that is not tabix'ed, you
				     must specify 'rsid_col_index', 'chr_col_index', and 'pos_col_index' in the config file")
			}

			if ("genome_build" %in% names(rsid_config))
			{
				print("Warning: you have specified both the 'genome_build' and the 'dbsnp_file' 
				      parameters for the 'add_rsids' tool. The 'genome_build' will be ignored for this
				      step and the custom rsid file used instead, regardless of its genome version.'")
			}
		} else 
		{
			
			rsid_config$use_tabix = "TRUE"
			rsid_config$rsid_col_index = 3
			if (!("genome_build" %in% names(rsid_config)))
			{
				stop("Config ERROR: must either specify an rsid mapping file, or
				     specify the 'genome_build' parameter")	
			}
			else if (!(rsid_config$genome_build %in% c("hg19", "hg38")))
			{
				stop("Config ERROR: 'genome_build' parameter for 'add_rsids' must equal 
				     either 'hg19' or 'hg38'")	
			}
			
			if (rsid_config$genome_build == "hg19")
			{
				rsid_config$dbsnp_file = "data/dbsnp/hg19/common_all_20170710.vcf.gz"
			} else if (rsid_config$genome_build == "hg38")
			{
				rsid_config$dbsnp_file = "data/dbsnp/hg38/common_all_20170710.vcf.gz"
			}

		}
		
		# If no output file specified, use a default one
		if (!("out_file" %in% names(rsid_config)))
		{
			rsid_config$out_file = "add_rsids_out.txt"
		}

		new_config$tool_settings$add_rsids = rsid_config
	}

	return(new_config)
}


