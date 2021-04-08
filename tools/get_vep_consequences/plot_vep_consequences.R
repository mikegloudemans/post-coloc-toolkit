require(rjson)
require(dplyr)
require(ggplot2)
require(reshape2)
require(cowplot)

plot_vep_consequences = function(config_file, coloc_file, annotations_file, output_directory)
{
	config = fromJSON(file=config_file)$plot_vep_consequences
	config$coloc_file = coloc_file
	config$annotations_file = annotations_file
	config$output_directory = output_directory

	
	# Load results, errors, skips files
	coloc_results = load_coloc_file(config)
	annotations_results = load_vep_file(config)
	
	# Filter down to colocalized loci...
	coloc_results = coloc_results %>% filter(coloc_status == "colocalized")
	coloc_veps = annotations_results %>% filter(locus %in% unique(coloc_results$locus))

	# How many loci had each annotation?
	annotation_counts = coloc_veps %>% group_by(consequence) %>% summarize(count = length(unique(locus)), percent = count / length(unique(coloc_veps$locus))*100)
	
	# Output as bar plot then, including one stratifying by consequence.
	
	annotation_counts$consequence = factor(annotation_counts$consequence,
				    levels = c("3PRIME_UTR", "5PRIME_UTR", "CANONICAL_SPLICE", "DOWNSTREAM", "INTERGENIC", "INTRONIC", 
					       "NONCODING_CHANGE", "NON_SYNONYMOUS", "REGULATORY", "SPLICE_SITE", "STOP_GAINED", 
					       "SYNONYMOUS", "UPSTREAM"),
				    labels = c("3' UTR", "5' UTR", "inside splice motif", "downstream", "intergenic", "intronic", "noncoding transcript", "non-synonymous", "regulatory element", "near splice motif", "stop-gain", "synonymous", "upstream"))

	annotation_counts$impact = ""
	annotation_counts$impact[annotation_counts$consequence %in% c("inside splice motif", "stop-gain")] = "high"
	annotation_counts$impact[annotation_counts$consequence %in% c("non-synonymous")] = "moderate"
	annotation_counts$impact[annotation_counts$consequence %in% c("near splice motif", "synonymous")] = "low"
	annotation_counts$impact[annotation_counts$consequence %in% c("3' UTR", "5' UTR", "downstream", "intergenic", "intronic", "noncoding transcript", "regulatory element","upstream")] = "modifier"
	annotation_counts$impact = factor(annotation_counts$impact, levels = c("high", "moderate", "low", "modifier"))

	impact_counts = data %>% group_by(impact) %>% summarize(count=length(unique(locus)), percent = count/length(unique(coloc_veps$locus))*100)
	impact_counts = data %>% filter(impact != "modifier")
	impact_counts = data %>% arrange(percent)

	plot_results(annotation_counts, impact_counts, config)

}

load_coloc_file = function(config)
{
	t = read.table(config$coloc_file, header=TRUE, sep="\t")

	return(t)	
}

load_vep_file = function(config)
{
	t = read.table(config$vep_file, header=TRUE, sep="\t")

	return(t)	
}

args = commandArgs(trailingOnly=TRUE)

config_file = args[1]
coloc_file = args[2]
annotations_file = args[3]
output_directory = args[4]
plot_vep_consequences(config_file, coloc_file, annotations_file, output_directory)

plot_results = function(effect_counts, impact_counts, config)
{
	g_d1 <- ggplot(data=effect_counts, aes(x=consequence, y=percent, fill=impact)) +
		geom_bar(stat="identity", color="black", position="dodge") +
		geom_text(color="black", aes(label=count, y=percent+6), size=3) +
		coord_flip() +
		theme_minimal() +
		scale_fill_manual(values = c("red", "yellow", "green", "grey85")) +
		theme(axis.title.y = element_blank())

	# Also include one with aggregate info
	g_d2_pre <- ggplot(data=impact_counts, aes(x=impact, y=percent, fill=impact)) +
		geom_bar(stat="identity", color="black", position="dodge") +
		geom_text(color="black", aes(label=count, y=percent+6), size=3) +
		coord_flip() +
		ylab("% of coloc loci with annotation") +
		ylim(c(0,100)) +
		theme_minimal() +
		scale_fill_manual(values = c("red", "yellow", "green", "grey85")) +
		theme(legend.position = "none") +
		theme(axis.title.y = element_blank())

	g_d2 = plot_grid(NULL, g_d2_pre, NULL, ncol=3, rel_widths=c(0.13,0.69,0.18))

	g_d = plot_grid(g_d1, g_d2, nrow=2, rel_heights=c(2.2,1))

	ggsave(config$output_file, height=4, width=15)
}
