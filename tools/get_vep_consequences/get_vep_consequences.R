require(reshape2)
require(ggplot2)
require(dplyr)
require(readr)
require(rjson)

# For each locus, make a list of possible causal variants (based on LD with lead variants or CLPP-mod threshold or something)
# and for each, list the variant's predicted effect from VEP

# CADD file, which contains VEP consequence annotations, the field we're looking for
cadd_file = "data/cadd/whole_genome_SNVs_inclAnno.tsv.gz"
# LD file, which contains all SNP pairs with r^2 >= 0.8
ld_file = "data/ld/EUR_geno_hg38.txt.gz"

#####################################################
### Setup and loading config settings
#####################################################

config_file = commandArgs(trailingOnly=TRUE)[1]

# Load pre-specified config file
config = fromJSON(file=config_file)

colocs = read.table(paste0(config$out_dir, "/clpp_results_categorized_", config$analysis_date, ".txt"), sep="\t", header=TRUE)

#####################################################
### Get list of variants to annotate
#####################################################

# Figure out which columns of the LD file we need
ld_head = unlist(strsplit(unlist(system(sprintf("zcat %s | head -n 1", ld_file), intern=TRUE)), "\t"))
ld_r2 = which(ld_head == "ld_r2")
ld_chr = which(ld_head == "snp2_chrom_hg38")
ld_pos = which(ld_head == "snp2_pos_hg38")

# Figure out which columns of the CADD file we need, to get the VEP consequence field
head = system(sprintf("zcat %s | tail -n +2 | head -n 1", cadd_file), intern=TRUE)
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
# For each locus
for (l in 1:max(colocs$locus))
{
	print(l)
	# Get all lead SNPs for that locus (regardless of GWAS trait)
	locus_colocs = colocs %>% filter(locus == l) %>% select(chr, pos)
	
	locus_snps = locus_colocs %>% filter(!duplicated(locus_colocs))
	
	# Allocate memory for list of all SNPs in LD
	max_num_snps = 100000	# Couldn't possibly be more than this could it?
	all_snps = as.data.frame(list(chr=rep(0, max_num_snps), pos=rep(0, max_num_snps), ld=rep(0, max_num_snps)))

	snp_num = 0
	# For each lead SNP, get its LD buddies
	for (snp in 1:dim(locus_snps)[1])
	{	
		snp_num = snp_num + 1
		
		# Add the lead SNP to the list
		all_snps$chr[snp_num] = locus_snps$chr[snp]
		all_snps$pos[snp_num] = locus_snps$pos[snp]
		all_snps$ld[snp_num] = 1
			
		# Get all SNPs that are LD buddies of any of those lead SNPs (r2>0.8) and add them to a "set" for the locus
		ld = unlist(system(sprintf("tabix %s chr%s:%s-%s", ld_file, as.character(locus_snps$chr[snp]), 
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
			all_snps$ld[snp_num] = buddy[ld_r2]
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
		out = unlist(system(sprintf("tabix %s %s:%s-%s", cadd_file, as.character(all_snps$chr[snp]), 
					    as.character(all_snps$pos[snp]), as.character(all_snps$pos[snp])), intern=TRUE))
		for (o in 1:length(out))
		{
			loc_table_idx = loc_table_idx + 1
			row = sapply(out[o], function(x) {as.character(strsplit(x, "\t")[[1]])})
			vep = row[vep_index]
			ref = row[ref_index]
			alt = row[alt_index]
			gene = row[gene_index]
			feature = row[feature_index]

			all_loci$chr[loc_table_idx] = all_snps$chr[snp]
			all_loci$pos[loc_table_idx] = all_snps$pos[snp]
			all_loci$r2[loc_table_idx] = all_snps$ld[snp]
			all_loci$locus[loc_table_idx] = l
			all_loci$vep[loc_table_idx] = vep
			all_loci$ref[loc_table_idx] = ref
			all_loci$alt[loc_table_idx] = alt
			all_loci$gene[loc_table_idx] = gene
			all_loci$feature[loc_table_idx] = feature
		}
	}
}

all_loci = all_loci[1:loc_table_idx,]
all_loci$id = paste(all_loci$vep, all_loci$chr, all_loci$pos, all_loci$ref, all_loci$alt, all_loci$gene, all_loci$feature, all_loci$r2, sep=".")
vep_summary = all_loci %>% group_by(locus) %>% summarize(paste(id, collapse="|"))

# Write output files

# One row per variant
write.table(all_loci, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, file=paste0(config$out_dir, "/vep_summary_melted_", config$analysis_date, ".txt"))
# One row per locus
write.table(vep_summary, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, file=paste0(config$out_dir, "/vep_summary_", config$analysis_date, ".txt"))

