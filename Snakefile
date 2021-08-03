
rule all:
	input:
		expand("output/make_heatmaps/{study}_completion_indicator.tmp", study=config["name"])

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

rule assign_locus_numbers:
	input:
		"output/add_hgnc_names/{study}_colocalization_results.txt"
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

rule make_heatmaps:
	input:
		"output/classify_results/{study}_colocalization_results.txt"
	output:
		"output/make_heatmaps/{study}_completion_indicator.tmp"
	params:
		config = "config/{study}.config",
		out_base = "output/make_heatmaps/{study}"
	shell:
		"Rscript tools/make_heatmaps/make_heatmaps.R {params.config} {input} {params.out_base} {output}"


