require(rjson)
require(dplyr)
require(ggplot2)
require(reshape2)
require(cowplot)

plot_category_bars = function(config_file, in_base, out_base, completion_indicator)
{
	config = fromJSON(file=config_file)$plot_category_bars
	config$in_base = in_base
	config$out_base = out_base
	config$completion_indicator = completion_indicator

	strata = Sys.glob(sprintf("%s_class_summary_*.txt", config$in_base))

	for (stratum in strata)
	{
		counts = load_summary_file(stratum)

		plot_results(counts, config)
	}

	system(sprintf("touch %s", config$completion_indicator))
}

load_summary_file = function(in_file)
{
	d = read.table(in_file, header=TRUE)
	return(d)
}

plot_results = function(data, config)
{
	datacopy = data

	class = colnames(datacopy)[1]

	colnames(datacopy)[1] = "stratum"

	# Re-order categories if order specified in config file
	if (("orders" %in% names(config)) && (class %in% names(config$orders)))
	{
		datacopy = datacopy %>% filter(stratum %in% config$orders[[class]])
		datacopy$stratum = factor(datacopy$stratum, levels = rev(config$orders[[class]]))
	}

			
	if (("colors" %in% names(config)) && (class %in% names(config$colors)))
	{
		if (!(("orders" %in% names(config)) && (class %in% names(config$orders))))
		{	
			sprintf("ERROR: If specifying a custom color scheme for %s, you also need to specify an order
				for the categories.", class)
			stop()
		}
		
		config$colors[[class]] = config$colors[[class]][config$orders[[class]] %in% datacopy$stratum]

		plot <- ggplot(data=datacopy, aes(x=stratum, y=num_loci, fill=stratum)) + 
			scale_fill_manual(values=rev(config$colors[[class]]))
	} else
	{
		plot <- ggplot(data=datacopy, aes(x=stratum, y=num_loci))
	}

	plot = plot + geom_bar(stat="identity", color="black", position=position_dodge())+
		geom_text(aes(label = num_loci, y = num_loci + 0.07 * max(datacopy$num_loci)-1))+
		coord_flip() +
		ylab("# loci") +
		theme_minimal() +
		theme(legend.position = "none") +
		theme(axis.title.y = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
			axis.text = element_text(size = 12))


	plot

	ggsave(sprintf("%s_%s_barplot.pdf", out_base, class), height=4, width=6)
}

args = commandArgs(trailingOnly=TRUE)
config_file = args[1]
in_base = args[2]
out_base = args[3]
completion_indicator = args[4]
plot_category_bars(config_file, in_base, out_base, completion_indicator)

