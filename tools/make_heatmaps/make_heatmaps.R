### Code by Brunilda Balliu 12/21/2019
# Updated by Mike Gloudemans 6/18/2020 and beyond

suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(require(rjson)))

chunk_size = 100
row_height=0.21
col_width=0.2

# Load default color scheme
source("tools/make_heatmaps/color_scheme.R")

############################################################
### Create colocalization heatmaps
############################################################

make_heatmaps = function(config_file, input_file, output_directory, completion_indicator)
{
	# Load config file
	config = fromJSON(file=config_file)$make_heatmaps
	config$input_file = input_file
	config$output_directory = output_directory
	
	if ("rows_per_page" %in% names(config))
	{
		chunk_size = as.numeric(config$rows_per_page)
	}

	# Load results table
	coloc_res = get_coloc_results(config$input_file, config)
	coloc_res = coloc_res %>% arrange(-score)
	
	# Make an individual split for every stratification wanted.
       	# Specify this in the config file	
	for (strat in config$file_strata)
	{

		coloc_res_tmp = coloc_res
		coloc_res_tmp$split_column = ""
		if(!("split_factors" %in% names(strat)))
		{
			coloc_res_tmp$split_column = "no-split"
			print("Plotting without stratification")
		} else
		{
			print(sprintf("Stratifying by %s", strat$split_factors))
		}

		for (column in strat$split_factors)	
		{
			if (sum(coloc_res_tmp$split_column != "") == 0)
			{
				coloc_res_tmp$split_column = coloc_res_tmp[[column]]
			}
			else
			{
				coloc_res_tmp$split_column = paste(coloc_res_tmp$split_column, coloc_res_tmp[[column]], sep="-")
			}
		}
		coloc_res_tmp$split_column = factor(coloc_res_tmp$split_column)

		if ("gwas_blacklist" %in% names(strat))
		{
			for (bl in strat$gwas_blacklist)
			{
				coloc_res_tmp = coloc_res_tmp %>% filter(gwas_label != bl)
			}
			coloc_res_tmp$gwas_label = factor(x=coloc_res_tmp$gwas_label)
		}

		coloc_res_tmp = collapse_axis_factors(coloc_res_tmp, config)
		
		# Classify the coloc results to determine what color they'll be in the plot
		coloc_res_tmp$coloc_class = get_heatmap_classes(coloc_res_tmp, config)

		# Find locus-gene pairs with no coloc at all
		# Label rows as "blank" if they have no matches, so we can leave them out of plots
		coloc_res_tmp$blanks = ""
		for (y_factor in unique(coloc_res_tmp$y_factor))
		{
			locus_matches = coloc_res_tmp[coloc_res_tmp$y_factor == y_factor,]
			if (sum(locus_matches$coloc_class!="none") == 0)
			{
				coloc_res_tmp$blanks[(coloc_res_tmp$y_factor == y_factor)] = "blank"
			}
			else
			{
				coloc_res_tmp$blanks[(coloc_res_tmp$y_factor == y_factor)] = "okay"
			}	
		}

		# Remove rows (gene-locus pairs) with no colocs at all, if desired
		if (lower(strat$concise) == "true")
		{
			coloc_res_tmp = coloc_res_tmp %>% filter(blanks != "blank")
		}

		# If there are multiple SNPs in the same plot cell, then
		# just pick the one with the highest CLPP mod
		coloc_res_tmp = coloc_res_tmp %>% arrange(-score)
		coloc_res_tmp = coloc_res_tmp[!duplicated(coloc_res_tmp[,c("x_factor", "y_factor")]),]

		coloc_res_tmp = coloc_res_tmp %>% arrange(y_factor)


		if (lower(config$cluster) == "true")
		{
			# Binarize cells into colocalize or non-colocalized
			coloc_res_tmp$clust_stat = 0
			coloc_res_tmp$clust_stat[coloc_res_tmp$coloc_class == "none"] = 1

			dc = dcast(coloc_res_tmp, y_factor ~ x_factor, value.var = "clust_stat")
			grid = as.matrix(dc[,-1])
			rownames(grid) = dc[,1]
			grid[is.na(grid)] = 0
			
			dst = suppressWarnings(dist(dc, method = "binary"))
			dst[is.na(dst)] = 1
			h = hclust(dst)

			coloc_res_tmp$y_factor = factor(coloc_res_tmp$y_factor, levels = rownames(grid)[h$order])
		}

		### Plot coloc results
		plot_heatmap(coloc_res_tmp, strat, config)
	}

	system(sprintf("touch %s", completion_indicator))

}

######################################################

# Function to create tileplot of coloc results
plot_coloc_results_function=function(data){
  plot=ggplot(data =  data,
              mapping = aes(x=x_factor, y = y_factor)) + 
    geom_tile(mapping = aes(fill=coloc_class),color="black") + 
    scale_fill_manual(name="colocalization score", values = color_scheme ,drop=FALSE) + 
    geom_text(mapping = aes(label=cross, color=cross_col), size=4) +
    scale_colour_manual(values=c("black", "white")) +
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(angle=90, size=15, hjust = 1, vjust=.5), 
          axis.text.y = element_text(size=12), 
          legend.position = 'top', 
          legend.box = "vertical", 
          legend.text = element_text(size=15), 
          legend.title = element_text(size = 12)) +   
	  guides(color=FALSE) + 
    scale_x_discrete(drop=FALSE)
  return(plot)
}

get_coloc_results = function(coloc_file, config)
{
	# Read coloc results
	coloc_res=as.data.frame(fread(file = coloc_file, sep = '\t', header = T, check.names = F))
	
	# Identify the QTL type and tissue for each coloc test
	coloc_res$qtl_type = coloc_res[[config[["type_column"]]]]
	coloc_res$tissue = coloc_res[[config[["tissue_column"]]]]
	coloc_res$gwas_label = coloc_res[[config[["gwas_column"]]]]
	
	# If tissue order not specified, just do it alphabetically
	if (!("tissue_order" %in% names(config))){
		config$tissue_order = sort(coloc_res$tissue)
	}
	coloc_res$tissue = factor(x=coloc_res$tissue, levels=config$tissue_order)

	# If tissue order not specified, just do it alphabetically
	if (!("gwas_order" %in% names(config))){
		config$gwas_order = sort(coloc_res$gwas_label)
	}
	coloc_res$gwas_label = factor(x=coloc_res$gwas_label, levels=config$gwas_order)

	if (("label_individual_cells" %in% names(config)) && (config$label_individual_cells == "True"))
	{
		# IF annotating with crosses and other marks to denote which ones 
		# aren't trustworthy, do it here.
		coloc_res = label_individual_cells(coloc_res)
	} else if (("put_scores_in_cells" %in% names(config)) && (config$put_scores_in_cells) == "True")
	{
		coloc_res = put_scores_in_cells(coloc_res, config)
	} else
	{
		coloc_res$cross = ""
		coloc_res$cross_col = ""
	}

	coloc_res$cross = factor(coloc_res$cross)
	coloc_res$cross_col = factor(coloc_res$cross_col)

	return(coloc_res)
}

put_scores_in_cells = function(coloc_res, config)
{
	# Annotate results with coloc scores
	# Will require further customization for any scores that aren't percentage based like CLPP or H4PP

	coloc_res_tmp = coloc_res
	coloc_res_tmp$cross = as.character(round(coloc_res_tmp$score, 2) * 100)
	coloc_res_tmp$cross_col = "black"
	coloc_res_tmp$cross_col[coloc_res_tmp$coloc_status == "coloc"] = "white"

	return(coloc_res_tmp)
} 

plot_heatmap = function(coloc_res, strat, config)
{
	out_sub_folder = strat$out_dir
	dir.create(paste0(config$output_directory, "/", out_sub_folder), recursive = TRUE, showWarnings=FALSE)
	if ("constrain_split" %in% names(strat))
	{
		splits = strat$constrain_split
	} else
	{
		splits = unique(levels(coloc_res[["split_column"]]))
	}
	for (i in 1:length(splits))
	{
		split_col = splits[i]
		
		print(paste0("Plotting by ", split_col))

		tmp_data=coloc_res %>% 
			filter(coloc_res[["split_column"]] == split_col)
		
		if (dim(tmp_data)[1] == 0)
		{
			next
		}

		row_count = length(unique(tmp_data$y_factor))

		max_chunk = floor(row_count / chunk_size + 1)

		for (chunk in 1:max_chunk)
		{
			tmp_chunk = tmp_data[tmp_data$y_factor %in% unique(tmp_data$y_factor)[(1+(chunk-1)*chunk_size):min(chunk_size*chunk, row_count)],]
			tmp_chunk$y_factor = factor(tmp_chunk$y_factor, levels = rev(unique(tmp_chunk$y_factor)))

			# Genes per locus
			genes_per_locus = tmp_chunk %>% group_by(locus) %>% summarize(genes_at_locus=length(unique(y_factor))) %>% arrange(locus)

			num_cols = length(levels(tmp_chunk$x_factor))
			if (("x_axis_collapse" %in% names(config)) && ((config$x_axis_collapse == "tissues") || (config$x_axis_collapse == "tissues-gwas")))
			{
				num_tissues = 1
			} else
			{
				num_tissues = length(levels(coloc_res$tissue))
			}
			num_vert_bars = num_cols / num_tissues - 1
			num_rows = length(unique(tmp_chunk$y_factor))

			if ((("cluster" %in% names(config)) && (lower(config$cluster) == "true")) ||
			    (("y_axis_collapse" %in% names(config)) && (config$y_axis_collapse == "genes")))
			{
				# It doesn't make sense to separate loci if they're clustered
				# It also doesn't make sense to separate loci if each row is a distinct locus
				num_horz_bars = 0
			} else 
			{
				num_horz_bars = length(unique(tmp_chunk$locus))
			}
			horz_breaks = cumsum(genes_per_locus$genes_at_locus)

			y_margin_approx_size = 0.2*max(c(nchar(as.character(unique(tmp_chunk$y_factor))), 0))
			x_margin_approx_size = 0.2*max(c(nchar(as.character(unique(tmp_chunk$x_factor))), 0))+1.3

			plot=plot_coloc_results_function(data = tmp_chunk)
			
			if (num_vert_bars != 0)
			{
				my.vertical.lines<-data.frame(x=seq(0, (num_vert_bars-1)*num_tissues, by = num_tissues) + num_tissues + 0.5, y = rep(0.5, num_vert_bars), 
					xend=seq(0, (num_vert_bars-1)*num_tissues, by = num_tissues) + num_tissues + 0.5, yend = rep(num_rows + 0.5, num_vert_bars))
				plot = plot + geom_segment(data=my.vertical.lines, aes(x,y,xend=xend, yend=yend), size=1, inherit.aes=F)
			}
			if (num_horz_bars != 0)
			{
				my.horizontal.lines<-data.frame(x=rep(0.5, num_horz_bars), y=num_rows-horz_breaks+0.5, 
					xend=rep(num_cols+0.5, num_horz_bars), yend=num_rows - horz_breaks+0.5)
				plot = plot + geom_segment(data=my.horizontal.lines, aes(x,y,xend=xend, yend=yend), size=0.25, inherit.aes=F)
			}
			plot

			ggsave(filename = paste0(config$output_directory, '/', out_sub_folder, '/CLPP_group_',split_col,'.part', chunk, '.pdf'), plot = plot, width = y_margin_approx_size+(col_width*num_cols), height = x_margin_approx_size+(row_height*num_rows), limitsize = F)
		}
	}
}

get_heatmap_classes = function(coloc_res, config)
{
        get_coloc_type = function(qtl_type, coloc_status, x_factor, y_factor)
	{
		qtl_type = tolower(qtl_type)
		#print(qtl_type)
		#print(sum(qtl_type == "eqtl"))
	      if ((sum((qtl_type == "eqtl") & (coloc_status == "coloc")) > 0) &&
		      (sum((qtl_type == "sqtl") & (coloc_status == "coloc")) > 0))
	      {
		      return("both")
	      }
	      if ((sum((qtl_type == "eqtl") & (coloc_status == "coloc")) > 0) &&
		      (sum((qtl_type == "sqtl") & (coloc_status == "coloc")) == 0))
	      {
		      return("eqtl")
	      }
	      if ((sum((qtl_type == "eqtl") & (coloc_status == "coloc")) == 0) &&
		      (sum((qtl_type == "sqtl") & (coloc_status == "coloc")) > 0))
	      {
		      return("sqtl")
	      }
	      return("none")
        }
	coloc_types = coloc_res %>% group_by(x_factor, y_factor) %>% summarize(coloc_class = get_coloc_type(qtl_type, coloc_status, x_factor, y_factor))

	coloc_res = full_join(coloc_res, coloc_types, by = c("x_factor", "y_factor"))

	coloc_res$coloc_class=factor(x = coloc_res$coloc_class, 
	levels=c("none", 
		 "sqtl",  
		 "eqtl",
		 "both"))
}

collapse_axis_factors = function(coloc_res_tmp, config)
{
	coloc_res = coloc_res_tmp

	if ("locus_selection_list" %in% names(config))
	{
		coloc_res = coloc_res %>% filter(locus %in% config$locus_selection_list)
	}
	if ("gene_selection_list" %in% names(config))
	{
		coloc_res = coloc_res %>% filter(ensembl %in% config$gene_selection_list)
	}

	# Collapse genes by locus
	# Remove all but the best coloc at each locus, regardless of gene
	if (("y_axis_collapse" %in% names(config)) && (config$y_axis_collapse == "genes"))
	{
		coloc_res$y_factor = coloc_res$locus
	} else
	{
		coloc_res$y_factor = paste(coloc_res$locus, coloc_res$hgnc, sep="-")
	}

	coloc_res = coloc_res %>% arrange(as.numeric(coloc_res$locus))
	coloc_res$y_factor = factor(coloc_res$y_factor, levels = unique(coloc_res$y_factor))
	
	# Collapse across tissues
	# Remove all but the best tissue for a given
	if (("x_axis_collapse" %in% names(config)) && (config$x_axis_collapse == "tissues"))
	{
		coloc_res$x_factor = coloc_res$gwas_label
	} else if (("x_axis_collapse" %in% names(config)) && (config$x_axis_collapse == "gwas"))
	{
		coloc_res$x_factor = coloc_res$tissue
	} else if (("x_axis_collapse" %in% names(config)) && (config$x_axis_collapse == "tissues-gwas"))
	{
		coloc_res$x_factor = "any"
		coloc_res$x_factor = factor(coloc_res$x_factor)
	} else
	{
		coloc_res$x_factor = paste(coloc_res$gwas_label, coloc_res$tissue, sep="-")
		coloc_res$x_factor = factor(coloc_res$x_factor, levels = as.vector(t(outer(levels(coloc_res$gwas_label), levels(coloc_res$tissue), FUN = "paste", sep="-"))))
	}
	# No need to sort on the x-axis, because the desired orders have already been pre-specified
	
	return(coloc_res)
}

args = commandArgs(trailingOnly=TRUE)

config_file = args[1]
input_file = args[2]
output_directory = args[3]
completion_indicator = args[4]

make_heatmaps(config_file, input_file, output_directory, completion_indicator)
