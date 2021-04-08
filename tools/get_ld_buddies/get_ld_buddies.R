require(rjson)
require(dplyr)

# For each locus, make a list of possible causal variants (based on LD with lead variants or CLPP-mod threshold or something)
# and for each, list the variant's predicted effect from VEP

get_ld_buddies = function(config_file, input_file, output_file)
{
	# Load config file
	config = fromJSON(file=config_file)$get_ld_buddies
	config$input_file = input_file
	config$output_file = output_file

	# Load results table
	results = load_results_file(config)

	# Get info on the LD file
	ld_info = get_ld_file_info(config)
	ld_r2 = ld_info$ld_r2
	ld_chr = ld_info$ld_chr
	ld_pos = ld_info$ld_pos

	# Initialize a data frame to hold all relevant SNPs
	max_total_snps = 1000000
	all_loci = as.data.frame(list(chr=rep("", max_total_snps), pos=rep(0, max_total_snps), locus=rep(0, max_total_snps), r2=rep("", max_total_snps)), stringsAsFactors=FALSE)

	loc_table_idx = 0
	# For each locus
	for (l in sort(unique(results$locus)))
	{
		print(l)

		# Get all colocalized lead SNPs for that locus (regardless of GWAS trait)
		locus_colocs = results %>% filter(locus == l) %>% select(chr, pos)
		locus_snps = locus_colocs %>% filter(!duplicated(locus_colocs))
		
		# Allocate memory for list of all SNPs in LD
		max_num_snps = 100000
		all_snps = as.data.frame(list(chr=rep(0, max_num_snps), pos=rep(0, max_num_snps), ld=rep(0, max_num_snps)))

		snp_num = 0
		# For each lead SNP within this locus, get its LD buddies
		for (snp in 1:dim(locus_snps)[1])
		{	
			snp_num = snp_num + 1

			# Add the lead SNP to the list
			all_snps$chr[snp_num] = locus_snps$chr[snp]
			all_snps$pos[snp_num] = locus_snps$pos[snp]
			all_snps$ld[snp_num] = 1
				
			# Get all SNPs that are LD buddies of any of those lead SNPs (r2>0.8) and add them to a "set" for the locus
			ld = unlist(system(sprintf("tabix %s chr%s:%s-%s", config$ld_file, as.character(locus_snps$chr[snp]), 
						   as.character(locus_snps$pos[snp]), as.character(locus_snps$pos[snp])), intern=TRUE))
			if (length(ld)==0)
			{
				next
			}
			for (i in 1:length(ld))
			{
				snp_num = snp_num + 1
				buddy = unlist(strsplit(ld[i], "\t"))
				all_snps$chr[snp_num] = gsub("chr", "", buddy[ld_chr])
				all_snps$pos[snp_num] = buddy[ld_pos]
				all_snps$ld[snp_num] = as.numeric(buddy[ld_r2])
			}
		}

		# Remove duplicates, since some SNPs might have been added more than once.
		# Keep the one with the highest LD with a lead GWAS SNP
		all_snps = all_snps[1:snp_num,]
		all_snps = all_snps %>% arrange(desc(ld))
		all_snps = all_snps %>% filter(!duplicated(paste(chr, pos)))

		# Now that we have all the SNPs and their LD buddies,
		# get VEP consequence annotations for the final table
		for (snp in 1:dim(all_snps)[1])
		{
			loc_table_idx = loc_table_idx + 1

			all_loci$chr[loc_table_idx] = all_snps$chr[snp]
			all_loci$pos[loc_table_idx] = all_snps$pos[snp]
			all_loci$r2[loc_table_idx] = all_snps$ld[snp]
			all_loci$locus[loc_table_idx] = l
		}
	}
	all_loci = all_loci %>% filter(!grepl("alt", chr) & !grepl("Un", chr) & !grepl("hap", chr)) %>% filter(chr != "")

	# One row per variant
	write.table(all_loci, file = config$output_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

get_ld_file_info = function(config)
{
	info = list()

	# Figure out which columns of the LD file we need
	ld_head = unlist(strsplit(unlist(system(sprintf("zcat %s | head -n 1", config$ld_file), intern=TRUE)), "\t"))
	info[["ld_r2"]] = which(ld_head == "ld_r2")
	if (config$genome == "hg38")
	{
		info[["ld_chr"]] = which(ld_head == "snp2_chrom_hg38")
		info[["ld_pos"]] = which(ld_head == "snp2_pos_hg38")
	} else if (config$genome == "hg19")
	{
		info[["ld_chr"]] = which(ld_head == "snp2_chrom_hg19")
		info[["ld_pos"]] = which(ld_head == "snp2_pos_hg19")
	} else
	{
		stop("ERROR: Please provide a valid setting for the 'genome' parameter
		     	in the config file.")
	}

	return(info)
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

get_ld_buddies(config_file, input_file, output_file)


