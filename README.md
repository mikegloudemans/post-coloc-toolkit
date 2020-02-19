# post-coloc-toolkit

Created by Mike Gloudemans

## Summary

This repository contains a suite of tools for exploring colocalization results.

The tools were designed to work with results in the output format of the pipeline
https://github.com/mikegloudemans/ensemble_colocalization_pipeline, but there's
no reason they won't work with the results of any quantitative colocalization pipeline 
if converted to the proper format.

Specific variations of these tools have been used in several collaborations; the aim 
of this package, however, is to enable generalized use of these tools with results
from _any_ GWAS and/or QTL study, with all key parameters optionally overridden by the
user. So, if you find ways in which the pipeline is incompatible with your own pipeline,
please feel free to open an issue or submit a pull request, and I'll do my best to
integrate your recommendations.

## Getting started

All modules in this package require a config file, which is specified as a command-line
parameter. Generally, it is simplest to share one single config file across all analyses.
Details on how to create such a file are given below in a separate section.

I've designed the tools so they can be run independently, but most often you'll want to
run them all together through a master script, which is implemented and described in a script in
the top-level folder, called `run_full_workflow.py`. Many of the steps in the workflow can
be skipped over by specifying this in the config file, which is described in the README
accompanying each individual step.


## Individual tools

Since the tools largely function independently of one another with just a few exceptions,
I've created an individual README for each one of the tools within its corresponding folder
describing its usage in detail. Here I provide a short description of what each of the
tools can do, and delineate how a typical workflow might look.

### cat-results

Concatenate individual results tables from single coloc runs into one big table for further
processing.

### post-hoc-filter

Optionally, filter out any studies or SNPs that meet exclusion criteria for this particular
post-colocalization analysis.

(e.g. Maybe you ran 10 loosely related GWAS traits, but you want to limit this analysis to 
type 2 diabetes loci with lead GWAS p-value < 1e-8, even though you performed colocalizations 
with a less stringent p-value threshold of p < 1e-5.)

### add-hgnc-names

Add columns to the results table showing HGNC gene names, when known. Not required for other analyses to run.

### add-rsids

Add columns to the results table showing HGNC gene names, when known. Not required for other analyses to run.

### assign-locus-numbers

Group nearby SNPs into clusters, and assign each group a locus number.

### classify-results

Sort candidate colocalization loci into individual categories

### last-qc-check

Verify that everything ran smoothly and that nothing unwanted made it through the filters.

### summarize-filter-counts

Summarize the number of SNPs, loci, and genes removed at each level of
filtering in the colocalization pipeline.

### plot-pipeline-filters

Generate a graphic showing the numbers of SNPs remaining after various 
pre- and post-colocalization filters.

Optionally, include additional final filters that have been specified
by the user through the classify-results tool.

### get-ld-buddies

### get-vep-consequences



## Config file

The config file is text file that represents a JSON object. (link to explanation of this format TODO)

The top level contains parameters shared across multiple modules such as the output directory, [analysis date?]
and location of raw colocalization input data. Most fields will be specified at this top level;
however, for parameters specific to a single module, they will be specified at a lower level as
described in the module's individual README.

An example minimal config file for full analysis is shown below. Additional
optional parameters are described in the READMEs for the relevant tools.

```
{
	# Example JSON
}
```


## Contact

