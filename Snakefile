rule all:
	input:
		"none"

# rule validate_config:
# Figure out later how to validate it before it goes in...	

# This might eventually become an optional step, but
# for now we're just assuming it's mandatory

# Filter down to the relevant files

rule post_hoc_filter:
	input:
		"data/coloc_results/{study}_colocalization_results.txt"
	output:
		"output/post_hoc_filter/{study}_colocalization_results.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/post_hoc_filter/post_hoc_filter.R {params.config} {input} {output}"

rule mutate_columns:
	input:
		"output/post_hoc_filter/{study}_colocalization_results.txt"
	output:
		"output/mutate_columns/{study}_colocalization_results.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/mutate_columns/mutate_columns.R {params.config} {input} {output}"


# This rule should probably be optional too
rule add_hgnc_names:
	input:
		"output/mutate_columns/{study}_colocalization_results.txt"
	output:
		"output/add_hgnc_names/{study}_colocalization_results.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/add_hgnc_names/add_hgnc_names.R {params.config} {input} {output}"

# This rule should probably be optional too
rule add_rsids:
	input:
		"output/add_hgnc_names/{study}_colocalization_results.txt"
	output:
		"output/add_rsids/{study}_colocalization_results.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/add_rsids/add_rsids.R {params.config} {input} {output}"

rule assign_locus_numbers:
	input:
		"output/add_rsids/{study}_colocalization_results.txt"
	output:
		"output/assign_locus_numbers/{study}_colocalization_results.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/assign_locus_numbers/assign_locus_numbers.R {params.config} {input} {output}"

rule classify_results:
	input:
		"output/assign_locus_numbers/{study}_colocalization_results.txt"
	output:
		"output/classify_results/{study}_colocalization_results.txt", summary = "output/classify_results/{study}_class_summary_completion_indicator.tmp"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/classify_results/classify_results.R {params.config} {input} {output} {output.summary}"

rule get_ld_buddies:
	input:
		"output/classify_results/{study}_colocalization_results.txt"
	output:
		"output/get_ld_buddies/{study}_ld_buddies.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/get_ld_buddies/get_ld_buddies.R {params.config} {input} {output}"

rule get_vep_consequences:
	input:
		"output/get_ld_buddies/{study}_ld_buddies.txt"
	output:
		"output/get_vep_consequences/{study}_vep_consequences.txt"
	params:
		config = "config/{study}.config"
	shell:
		"python tools/get_vep_consequences/get_vep_consequences.py {params.config} {input} {output}"

rule make_heatmaps:
	input:
		"output/classify_results/{study}_colocalization_results.txt"
	output:
		"output/make_heatmaps/{study}"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/make_heatmaps/make_heatmaps.R {params.config} {input} {output}"

rule plot_category_bars:
	input:
		"output/classify_results/{study}_class_summary_completion_indicator.tmp"
	output:
		"output/plot_category_bars/{study}_completion_indicator.tmp"
	params:
		config = "config/{study}.config", 
		in_base = "output/classify_results/{study}", 
		out_base = "output/plot_category_bars/{study}"

	shell:
		"Rscript tools/plot_category_bars/plot_category_bars.R {params.config} {params.in_base} {params.out_base} {output}"

