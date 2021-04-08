require(reshape2)
require(ggplot2)
require(dplyr)
require(readr)
require(rjson)

# For each LD-linked variant with top GWAS SNPs, list the variant's predicted effect from VEP

# TODO: Do it in batches by locus ... it's way to slow right now!

get_vep_consequences = function(config_file, input_ld_file, output_file)
{
	# Load config file
	config = fromJSON(file=config_file)$get_vep_consequences
	config$input_ld_file = input_ld_file
	config$output_file = output_file

	ld = load_ld_file(config)


	# Figure out which columns of the CADD file we need, to get the VEP consequence field
	skip_lines = 1
	if ("skip_cadd_lines" %in% names(config))
	{
		skip_lines = config$skip_cadd_lines
	}

	# TODO: Make this part more generalized to other consequence files
	head = system(sprintf("zcat %s | tail -n +%d | head -n 1", config$cadd_file, skip_lines+1), intern=TRUE)
	head = unlist(strsplit(head, "\t"))

	vep_index = which(head == "Consequence")
	ref_index = which(head == "Ref")
	alt_index = which(head == "Alt")
	gene_index = which(head == "GeneID")
	feature_index = which(head == "FeatureID")	

	# Initialize a data frame to hold all relevant SNPs
	max_total_snps = 1000000
	all_loci = as.data.frame(list(chr=rep(0, max_total_snps), pos=rep(0, max_total_snps), locus=rep(0, max_total_snps), vep = rep("", max_total_snps), ref=rep("", max_total_snps), alt=rep("", max_total_snps), gene=rep("", max_total_snps), feature=rep("", max_total_snps), r2=rep("", max_total_snps)), stringsAsFactors=FALSE)

	loc_table_idx = 0
	#for (i in 1:dim(ld)[1])
	for (i in 1:10)
	{
		print(i)
		snp = ld[i,]

		out = system(sprintf("tabix %s %s:%s-%s", config$cadd_file, as.character(snp$chr), 
					    as.character(snp$pos), as.character(snp$pos)), intern=TRUE)
	
		for (o in 1:length(out))
		{
			loc_table_idx = loc_table_idx + 1
			row = sapply(out[o], function(x) {as.character(strsplit(x, "\t")[[1]])})
			vep = row[vep_index]
			ref = row[ref_index]
			alt = row[alt_index]
			gene = row[gene_index]
			feature = row[feature_index]

			all_loci$chr[loc_table_idx] = snp$chr
			all_loci$pos[loc_table_idx] = snp$pos
			all_loci$r2[loc_table_idx] = snp$r2
			all_loci$locus[loc_table_idx] = snp$locus
			all_loci$vep[loc_table_idx] = vep
			all_loci$ref[loc_table_idx] = ref
			all_loci$alt[loc_table_idx] = alt
			all_loci$gene[loc_table_idx] = gene
			all_loci$feature[loc_table_idx] = feature
		}
	}

	all_loci = all_loci[1:loc_table_idx,]
	print(head(all_loci))
	print(tail(all_loci))

	# Write output files
	# One row per variant
	write.table(all_loci, file = config$output_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)       
}

load_ld_file = function(config)
{
	ld = read.table(config$input_ld_file, header=TRUE)

	return(ld)
}

args = commandArgs(trailingOnly=TRUE)

config_file = args[1]
input_ld_file = args[2]
output_file = args[3]

get_vep_consequences(config_file, input_ld_file, output_file)

