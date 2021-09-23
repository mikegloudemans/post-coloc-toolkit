suppressWarnings(suppressMessages(require(rjson)))
suppressWarnings(suppressMessages(require(dplyr)))

############################################################
### Add HGNC gene names to Ensembl-annotated file
############################################################

add_hgnc_names = function(config_file, input_file, output_file)
{
	config = fromJSON(file=config_file)$add_hgnc_names
	config$input_file = input_file
	config$output_file = output_file

	# Load results table
	results = load_results_file_for_hgnc(config)

	# Get table of mappings to HGNC genes
	genes = get_hgnc_table(config)

	# Remove entries with duplicated Ensembl IDs
	genes = genes[!duplicated(genes$ensembl),]
	
	# Combine HGNC table with results table
	results = merge(genes, results, by="ensembl", all.y=TRUE)
	results$hgnc = as.character(results$hgnc)
	results$ensembl = as.character(results$ensembl)
	results$hgnc[is.na(results$hgnc)] = results$ensembl[is.na(results$hgnc)]
	results$hgnc[results$hgnc == ""] = results$ensembl[results$hgnc == ""]

	write.table(results, config$output_file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
}

get_hgnc_table = function(config)
{

	# Load the table of mappings from Ensembl to HGNC	
	column_indices = as.numeric(c(config$ensembl_col_index, config$hgnc_col_index))
	genes = read.table(config$ensembl_to_hgnc_map_file, header=TRUE, sep="\t")

	genes = genes[,column_indices]
	colnames(genes) = c("ensembl", "hgnc")
	return(genes)
}

load_results_file_for_hgnc = function(config)
{
	data = read.table(file=config$input_file, header=TRUE)

	# Quick input check
	if (!("ensembl" %in% colnames(data)))
	{
		stop("input error: the input coloc results table must have an 'ensembl' column showing Ensembl IDs")
	}

	data$ensembl = substring(data$ensembl, 1, 15)

	return(data)
}

args = commandArgs(trailingOnly=TRUE)

config_file = args[1]
input_file = args[2]
output_file = args[3]
add_hgnc_names(config_file, input_file, output_file)
