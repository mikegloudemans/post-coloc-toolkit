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
		"output/classify_results/{study}_colocalization_results.txt", summary = "output/classify_results/{study}_class_summary.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/classify_results/classify_results.R {params.config} {input} {output} {output.summary}"

rule get_ld_buddies:
	input:
		"output/classify_results/{study}_colocalization_results.txt"
	output:
		"output/assign_locus_numbers/{study}_ld_buddies.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/get_ld_buddies/get_ld_buddies.R {params.config} {input} {output}"

