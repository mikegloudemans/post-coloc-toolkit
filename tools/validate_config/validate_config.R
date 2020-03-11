validate_config = function(config)
{ 
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
		if (!("raw_coloc_output_dirs" %in% names(cat_config)))
		{
			stop("Config ERROR: You must specify 'raw_coloc_output_dirs' in the 'cat_results' object, 
			     or else skip the cat_results step.")
		}
	}
}


