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

# This rule should probably be optional too
rule add_hgnc_names:
	input:
		"output/post_hoc_filter/{study}_colocalization_results.txt"
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

# This rule should probably be optional too
rule assign_locus_numbers:
	input:
		"output/add_rsids/{study}_colocalization_results.txt"
	output:
		"output/assign_locus_numbers/{study}_colocalization_results.txt"
	params:
		config = "config/{study}.config"
	shell:
		"Rscript tools/assign_locus_numbers/assign_locus_numbers.R {params.config} {input} {output}"

