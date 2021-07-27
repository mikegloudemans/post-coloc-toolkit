# post-coloc-toolkit

Created by Mike Gloudemans

## Summary

This repository contains a suite of tools for exploring colocalization results.

The tools were designed to work with results in the output format of the pipeline
https://github.com/mikegloudemans/ensemble_colocalization_pipeline, but
they will also work with the results of any quantitative gene prioritization pipeline 
if the appropriate column headers are included.

The aim of this package, is to enable generalized use of these tools with results
from _any_ GWAS and/or QTL study, with all key parameters optionally overridden by the
user. So, if you find ways in which the pipeline is incompatible with your own pipeline,
please feel free to a pull request, so that we might integrate your changes.

## Installation and demo

This pipeline is implemented using `snakemake`. 
TODO TODO how to install

How to run the demo

## Getting started

The settings for this toolkit are defined within a single config file, which is specified as a command-line
parameter. Details on how to create such a file are given below in a separate section.

The individual tools can be run independently, but most often you'll want to
run them all together, which is the default snakemake workflow if no other options
are specified.

## Typical workflow

The tools largely function independently of one another, with just a few dependencies.
Here I briefly describe each tool and delineate how a typical workflow might look. I then describe each rule
in greater detail, with parameter definitions and instructions for creating a custom config file.
These tools are listed in the order they will run for a typical workflow.

### post-hoc-filter

Filter out any studies or SNPs that meet exclusion criteria for this particular
post-colocalization analysis. In some cases, there might also be no required filters.

(e.g. Maybe you ran colocalization with adipose eQTLs and 10 loosely related GWAS traits, 
but you want to limit this analysis to only type 2 diabetes loci with lead eQTL p-value < 1e-10, 
even though you performed colocalizations with a less stringent p-value threshold of p < 1e-5.)

### mutate-columns

Change the names of columns, or add new columns based on a custom mapping from 
old column names. An example use might be specifying which 
QTL files correspond to effects in a certain tissue, or to eQTLs vs. sQTLs or other QTL types.

### add-hgnc-names

Add columns to the results table showing HGNC gene names, when known.

### assign-locus-numbers

Group nearby SNPs into discrete loci, and assign each group a locus number.

### classify-results

Organize loci by number of candidate genes and colocalizations, tissue(s) of colocalization,
trait(s) of colocalization, or other categorical filters to aid in gene prioritization.

### make-heatmaps

Make heatmaps in the style of (paper TBD...), illustrating colocalization results across
multiple dimensions. Results can be sorted and displayed using the categories determined
in the `classify-results` step.

## Creating a config file

The config file is simply a text file that represents a [JSON object](https://en.wikipedia.org/wiki/JSON#Syntax).
A JSON object is structurally (and syntactically) similar to a nested set of Python dictionaries and lists, and
I'll describe them as such here for simplicity.

An example config file is given at `config/examples/ir.config`. As seen in this file, the config file
is in essence a collection of dictionaries, one for each tool in the pipeline, specifying the parameters for running
that tool. (A further description of these parameters for each method is given later in the "Individual tool parameters"
section of the README.)

Example:

```
{
	"post_hoc_filter":
	{
		... (parameters for post_hoc_filter tool)
	},
	"mutate_columns":
	{
		... (parameters for mutate_columns tool)
	},
	"add_hgnc_names":
	{
		... (parameters for add_hgnc_names tool)
	},
	... (parameters for other tools)
}
```

As is hopefully evident in the above example, the config object is really nothing more than a 
nested collection of mini-config objects, one for each tool being run.

However, in case this still seems daunting, we have provided a script "wizard.py" that you
run to aid in creating your own config file. To run this script, execute in the command line:

```
python wizard.py
```

and follow the directions given.

^ NOTE: I have not yet completed this yet, so don't run `wizard.py` yet. Sorry.

## Input file format

The starting file is a TSV-formatted (tab-separated values) text file in which each
row represents the results of a single colocalization test. The name of this file should
correspond to the name of the associated config file; specifically, the config file
`config/{descriptor}.config` should be accompanied by an input data file called
`data/coloc_results/{descriptor}_colocalization_results.txt`

E.g. if the config file is `config/ir.config`,
the input data file should be located at `data/coloc_results/ir_colocalization_results.txt`.

The following columns are required in the file, and must be denoted in a header line.

* `ref_snp`: The central SNP used as a seed for selecting a colocalization window, usually the lead GWAS SNP.
* `gwas_trait`: The GWAS file being tested.
* `qtl_file`: The QTL file being tested.
* `feature`: The feature being tested in the QTL file; a gene for eQTLs, a splice junction for sQTLs, or something else for other QTL types.
* `n_snps`: The number of SNPs at the tested locus. Only required if results will be filtered based on this column.
* `score`: A numerical score representing the likelihood of colocalization. The scale is unimportant, but must be monotonic increasing (greater score means greater colocalization probability) and often includes values ranging from 0 to 1.
* `neg_log_gwas_pval`: The maximum (negative log p-value) of any GWAS SNP tested the locus. Only required if results will be filtered based on this column.
* `neg_log_qtl_pval`: The maximum (negative log p-value) of any QTL SNP tested the locus. Only required if results will be filtered based on this
column.
* `ensembl`: May be the same as `feature`, but this column maps the feature to a specific Ensembl gene.

## Included input and data files

A few subfolders of the `data` directory containing small text files are included as part of the repository
(you don't need to do anything to use them by default):

* `data/hgnc/ensembl_to_hgnc.txt` is a mapping of Ensembl IDs to HGNC gene IDs. You can substitute this file
	with your own mapping if you'd like; see the section on `add_hgnc_names` to learn more.
* `data/ldetect` contains mappings of genomic intervals to consistent locus numbers, with one file for the hg19 build
	and another for hg38. These intervals were determined with LDetect in EUR populations (Berisa & Pickrell 2016).
	Custom locus definitions or loci derived from other populations can be substituted here; however, all genomic
	positions must fall within one and only one locus for a mapping to be valid.
* `data/coloc_results` contains the results from a few colocalization runs used as examples and for the demos.

## Individual tool parameters

### `post_hoc_filter`

The post_hoc_filter tool is used to apply filters that were not done before colocalization testing.
Such filters might include
* A list of GWAS or QTL files to be excluded (or included) in subsequent analysis
* Maximum P-value thresholds for GWAS or QTL hits
* Minimum number of SNPs to be tested at a locus

The product of this script is to produce a filtered version of the input colocalization file,
but one in which colocalization tests not passing the filters are removed.

NOTE: This step is _not_ intended to remove all tests failing to pass a colocalization threshold.
It is important to hold on to such results for now, since they represent negative results. However,
we do label tests as colocalized or non-colocalized here based on a desired threshold, but without
performing any filtering based on these.

All parameters:

#### `colocalization_threshold`
A value indicating what level of the `score` column will indicate a colocalization when exceeded.

#### `kept_gwas`
Optionally, a list of values of the `gwas_trait` column that should be retained after filtering.
(All other values will be removed in the output.)

#### `removed_gwas`
Optionally, a list of values of the `gwas_trait` column that should be removed in filtering.
(All other values will be retained in the output.)

#### `kept_qtl`
Optionally, a list of values of the `qtl_file` column that should be retained after filtering.
(All other values will be removed in the output.)

#### `removed_qtl`
Optionally, a list of values of the `qtl_file` column that should be removed in filtering.
(All other values will be retained in the output.)

#### `gwas_pval_threshold`
An object indicating a "standard" maximum GWAS p-value to retrain and, optionally, "exceptions" of
alternate p-value thresholds for other traits (this latter option may be used if certain traits
have disproportionately high or low power to detect significant hits, but should
be always be used with caution and skepticism when relaxing threshold stringency). 
See the example config file below for a demonstration of how to use these values.

#### `qtl_pval_threshold`

Same as `gwas_pval_threshold`, just for QTL p-value instead of GWAS p-value.

#### `snp_count_min_threshold`

A minimum value to allow in the `n_snps` column, used for filtering colocalization tests
in which they're weren't enough SNPs profiled to make a confident prediction. The optimal
setting for this value depends on a lot of factors, but we generally recommend requiring at
least 20 SNPs to be tested per locus.

#### Example

An full example from the `ir.config` file:

```
{
	...,
        "post_hoc_filter":
        {
                "colocalization_threshold": 0.35,
                "kept_gwas":
                [
                        "BMI_GIANT_2018.txt.gz",
                        "T2D_Mahajan_Europeans.txt.gz",
                        "MI_adjBMI_European.txt.gz",
                        "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz",
                	...
		],
                "kept_qtl":
                [
                        "data/eqtls/gtex_v8/Adipose_Visceral_Omentum.allpairs.txt.gz.eQTLs.txt.gz",
                        "data/eqtls/gtex_v8/Adipose_Subcutaneous.allpairs.txt.gz.eQTLs.txt.gz",
                        "data/sqtls/gtex_v8/Adipose_Visceral_Omentum.sQTLs.txt.gz",
                        "data/sqtls/gtex_v8/Adipose_Subcutaneous.sQTLs.txt.gz",
			...
                ],
                "gwas_pval_threshold":
                {
                        "standard": 5e-8,
                        "exceptions":
                        {
                                "MI_adjBMI_European.txt.gz": 1e-5,
                                "MAGIC_ISI_Model_2_AgeSexBMI.txt.txt.gz": 1e-5
                        }
                },
                "qtl_pval_threshold":
                {
                        "standard": 1e-5
                },
                "snp_count_min_threshold": 20
        },
	...
}
```

### `mutate_columns`

The `mutate_columns` tool creates new columns in the file by mapping values of old columns to new ones.
The mapping can be many-to-one; that is, multiple 

Some examples of useful remappings might be to create columns that distinguish eQTL vs. sQTL,
tissue of QTL colocalization regardless of eQTL vs. sQTL, or higher-level trait categories such
as anthropometric traits vs. insulin/glucose measurments vs. others. Then these columns
can later be used for categorization of loci based on which types of colocalizations they contain
(e.g. eQTL vs. sQTL vs. both) or for sorting loci in an intuitive way when producing output heatmaps.

Alternatively, these columns can also be used simply to designate display names for files that
are better suited for showing as text in plots (e.g. "BMI_GIANT_2018.txt.gz": "BMI")

No existing columns will be dropped by this tool; only new columns will be created. The relevant parameter
is `mutations`, a list of dictionaries, each specifying an input (existing) column name, an output (non-existing)
column name, and the mapping of existing input entries to output entries.

#### Example

The following example from `config/ir.config` creates a column with shortened GWAS trait names for display,
one that specifies the tissue for each QTL file, and one that specifies the QTL type for each
QTL file. Note that multiple values in the input column can map to the same output column, which
may often be useful.

```
{
	...,
        "mutate_columns":
        {
                "mutations":
                [
                        {
                                "in": "gwas_trait",
                                "out": "gwas_short_trait",
                                "map":
                                {
                                        "BMI_GIANT_2018.txt.gz": "BMI",
                                        "FastGlu_MAGIC_Europeans.txt.gz": "FG",
                                        "FastInsu_adjBMI_MAGIC_Europeans.txt.gz": "FI",
                                        "T2D_Mahajan_Europeans.txt.gz": "T2D",
                                        "Type-2-Diabetes_Xue_2018.txt.gz": "T2D",
                                        "Type-2-Diabetes_Suzuki_2019.txt.gz": "T2D",
                                        "Type-2-Diabetes_Spracklen_2020.txt.gz": "T2D",
					...
                                }
                        },
                        {
                                "in": "qtl_file",
                                "out": "qtl_tissue",
                                "map":
                                {
                                        "data/eqtls/gtex_v8/Adipose_Visceral_Omentum.allpairs.txt.gz.eQTLs.txt.gz": "Visceral-Adipose",
                                        "data/eqtls/gtex_v8/Adipose_Subcutaneous.allpairs.txt.gz.eQTLs.txt.gz": "Subcutaneous-Adipose",
                                        "data/sqtls/gtex_v8/Adipose_Visceral_Omentum.sQTLs.txt.gz": "Visceral-Adipose",
                                        "data/sqtls/gtex_v8/Adipose_Subcutaneous.sQTLs.txt.gz": "Subcutaneous-Adipose",
					...
                                }
                        },
                        {
                                "in": "qtl_file",
                                "out": "qtl_type",
                                "map":
                                {
                                        "data/eqtls/gtex_v8/Adipose_Visceral_Omentum.allpairs.txt.gz.eQTLs.txt.gz": "eQTL",
                                        "data/eqtls/gtex_v8/Adipose_Subcutaneous.allpairs.txt.gz.eQTLs.txt.gz": "eQTL",
                                        "data/sqtls/gtex_v8/Adipose_Visceral_Omentum.sQTLs.txt.gz": "sQTL",
                                        "data/sqtls/gtex_v8/Adipose_Subcutaneous.sQTLs.txt.gz": "sQTL",
					...
                                }
                        }
                ]
        },
	...
}
```


### `add_hgnc_names`

`add_hgnc_names` takes the `ensembl` column from the input file and maps it to an "hgnc" column with
the corresponding standard HGNC gene name, when one exists.

The following parameters are required:

#### `ensembl_to_hgnc_map_file`

A file containing the mapping from Ensembl gene names to HGNC gene names. 

A default file "data/hgnc/ensembl_to_hgnc.txt" is supplied as part of the 
Github repository; however, other mappings can be obtained through
BioMart or other sources and included here as a TSV file with headers.

#### `ensembl_col_index`

The numerical (1-based) index of the column in `ensembl_to_hgnc_map_file` that
contains the Ensembl gene names.

#### `hgnc_col_index`

The numerical (1-based) index of the column in `ensembl_to_hgnc_map_file` that
contains the HGNC gene names.

#### Example

```
{
	...,
        "add_hgnc_names":
        {
                "ensembl_to_hgnc_map_file": "data/hgnc/ensembl_to_hgnc.txt",
                "ensembl_col_index": 1,
                "hgnc_col_index": 3
        },
	...
}
```

### `assign_locus_numbers`

This tool groups chromosomally nearby colocalization tests, including tests on different features
and/or traits, into predefined loci determined in European reference populations during the LDetect
algorithm. Optionally, other partitions of loci determined from other populations may also be used for this step,
but this currently would require you to edit the source code or replace the source files in `data/ldetect/` directly.

The only required parameter for this step is:

#### `genome_build`

Set to `hg19` or `hg38` to indicate which reference build of the genome was used for colocalization input files.

#### Example

```
{
	...,
        "assign_locus_numbers":
        {
                "genome_file": "hg38"
        },
	...
}
```

### `classify_results`

This important step can be used to classify entire loci based on the circumstances under which
colocalization was observed. Any number of classification "rules" can be specified within a dictionary
"rules" that maps rule names to descriptions about how to classify loci using these rules.

Note that these rules are applied not to each individual colocalization test, but collectively to all
colocalization tests performed at a single locus.

Here's a quick example of how that looks at the high level.

#### Example

```
{
        "classify_results":
        {
                "rules" :
                {
                        "rule1":
                        {
				...
			},
			"rule2":
			{
				...
			},
			...
			"rule99":
			{
				...
			}
		}
	}

```

There are two general types of rules that can be defined.
`num_colocs` rules are based on the number of candidate genes with overlapping QTLs
before running colocalization, and how many of these remain after running colocalization tests.
`specificity` rules are more general, and allow classification of loci with colocalization only
in certain subsets of tissues, traits, QTL types, or any other variable of interest.

Each rule will add an additional column to the output file that indicates in which category 
each locus belongs. A locus cannot belong to more than one class, and an error will be
raised if the classes are defined in a way that overlaps. (If a locus seems like it should
belong to more than one class, consider instead making two separate rules for the locus.)
Some loci may not belong to any class, although we generally recommend assigning each locus
to exactly one class, even if just named "None" or "Other".

#### `num_colocs` rules

These rules sort loci based on the number of candidate genes before running colocalization
and the number of genes passing the colocalization threshold. A reasonable default way of
defining these classes is shown in the example file `config/ir.config`, but other custom 
class definitions are possible too.

This type of rule is specified with a "type":"num_colocs" parameter in the rule object.
Then, a "categories" object defines an arbitrary number of discrete categories, each of which
is itself an object that defines the selection criteria for the category.

The allowable criteria are the following:
* `num_candidates_equals`
* `num_colocs_equals`
* `num_candidates_greater_than`
* `num_colocs_greater_than`
* `num_candidates_less_than`
* `num_colocs_less_than`

Criteria are defined as parameters within the rule(s), and each applies the exact filter
that it sounds like it does. `num_colocs` in the above criteria refers to the number of 
unique genes with at least one colocalization test passing the preivously specified colocalization
threshold. `num_candidates` refers to the number of unique genes that have a QTL overlapping the
GWAS hit and therefore were tested for colocalization at all, regardless of whether colocalization
was actually observed. 

#### Example:

This is the first of 3 rules in the example file `config/ir.config`. It
defines 5 different mutually distinct classes, and together these classes
encompass all possible loci.

```
{
	"classify_results":	
        {
                "rules" :
                {
                        "step1_coloc_sorting":
                        {
                                "type": "num_colocs",
                                "categories":
                                {
                                        "one_candidate_one_coloc":
                                        {
                                                "num_candidates_equals": "1",
                                                "num_colocs_equals": "1"
                                        },
                                        "one_candidate_no_coloc":
                                        {
                                                "num_candidates_equals": "1",
                                                "num_colocs_equals": "0"
                                        },
                                        "multi_candidate_one_coloc":
                                        {
                                                "num_candidates_greater_than": "1",
                                                "num_colocs_equals": "1"
                                        },
                                        "multi_candidate_multi_coloc":
                                        {
                                                "num_candidates_greater_than": "1",
                                                "num_colocs_greater_than": "1"
                                        },
                                        "multi_candidate_no_coloc":
                                        {
                                                "num_candidates_greater_than": "1",
                                                "num_colocs_equals": "0"
                                        }
                                }
                        },
			...
		}
	}
}
```

#### `specificity` rules

The other type of rule is `specificity` rules, which offer a broader range of custom classes.
These rules are indicated using the parameter "type":"specificity" in the rule object, and
also by selecting a particular existing `column` on which the rule is applied. This column
can be any column already existing in the input file, but will often be the column specifying
the GWAS trait or category, the QTL tissue, or the QTL type (eQTL, sQTL, or others).

The `specificity` rules also define a `categories` field, but unlike the `num_colocs` rules,
this field points to a list / array rather than an object / dictionary. The exact order in which
categories are specified here is important, because the tool will execute each rule, one at a time,
on only the loci for which a category has not already been assigned. That is, a locus may match
multiple categories, and if so, only the one appearing first in the list will be assigned to that locus.

These rules operate only on colocalized loci, and assign them to categories based on the range of values
observed in the field specified by the `column` parameter for _colocalized results only_. (E.g. in the example
that is shown below, a locus that was tested for colocalization only in pancreas tissue, and had no colocalization
for this test, will not match the "Pancreas" category's criteria.
Each category must assign a `class_name` (the value that will appear in the output field for loci matching 
this category's criteria) and any number of additional criteria.

The allowable criteria for this type of rule are the following:
* `contains_exactly`: Passes filter if the `column` contains colocalizations for all of the fields in this list, and no others.
* `contains_all`: Passes filter if the `column` contains colocalizations for all of the fields in this list, and maybe others.
* `contains_some`: Passes filter if the `column` contains colocalizations for at least one of the fields in this list, but not necessarily all.
* `contains_none`: Passes filter if the `column` contains no colocalizations for any of the fields in this list.
* `contains_only`: Passes filter if the `column` contains no colocalizations for any fields except those in this list, but not necessarily any colocalization for all or even any of the fields in this list.

#### Example 1: Tissue specificity

In this example, we use the `contains_exactly` keyword to identify loci that only have colocalizations 
in one tissue and no others. All other loci having colocalizations but not matching a particular will be labeled "Shared",
the last rule. Note that all loci in other categories also match the criteria for the "Shared" category, but they
will be assigned to those other categories instead because the "Shared" category is listed last and therefore is of
the lowest priority.

There is also one additional category that includes loci with both subcutaneous _and_ visceral adipose colocalizations,
but no colocalizations in any other tissues.

```
{
        "classify_results":
        {
                "rules" :
                {
                        "step2_tissue_sorting":
                        {
                                "type": "specificity",
                                "column": "qtl_tissue",
                                "categories":
                                [
                                        {
                                                "class_name": "AdpS",
                                                "contains_exactly":
                                                [
                                                        "Subcutaneous-Adipose"
                                                ]
                                        },
                                        {
                                                "class_name": "AdpV",
                                                "contains_exactly":
                                                [
                                                        "Visceral-Adipose"
                                                ]
                                        },
                                        {
                                                "class_name": "AdpS+V",
                                                "contains_exactly":
                                                [
                                                        "Subcutaneous-Adipose",
                                                        "Visceral-Adipose"
                                                ]
                                        },
                                        {
                                                "class_name": "Musk",
                                                "contains_exactly":
                                                [
                                                        "Skeletal-Muscle"
                                                ]
                                        },
                                        {
                                                "class_name": "Liver",
                                                "contains_exactly":
                                                [
                                                        "Liver"
                                                ]
                                        },
                                        {
                                                "class_name": "Pancreas",
                                                "contains_exactly":
                                                [
                                                        "Pancreas"
                                                ]
                                        },
                                        {
                                                "class_name": "Shared"
                                        }
                                ]
                        }
		}
	}
}
```
		
#### Example 2: Prioritizing specific traits

In contrast with the previous example, in which none of the category definitions but the "Shared" category
overlap one another, in this example a single loci could in principle meet the criteria for all
the categories, having at least one colocalization in each of 5 different defined priority groups of
traits. However, the order in which the traits are listed means that the "MI" and "ISI" colocalizations
will have the highest priority, followed by "T2D", "FG", and "FI" colocalizations, and all the way
down to "BMI" colocalizations, which will be classified as such only if there are no other colocalizations
with higher-priority traits.

```
{
	"classify_results":
	{
		"rules":
		{
                        "step3_gwas_sorting":
                        {
                                "type": "specificity",
                                "column": "gwas_short_trait",
                                "categories":
                                [
                                        {
                                                "class_name": "IR",
                                                "contains_some":
                                                [
                                                        "MI",
                                                        "ISI"
                                                ]
                                        },
                                        {
                                                "class_name": "T2D-FG-FI",

                                                "contains_some":
                                                [
                                                        "T2D",
                                                        "FG",
                                                        "FI"
                                                ]
                                        },
                                        {
                                                "class_name": "WHR",
                                                "contains_some":
                                                [
                                                        "WHR"
                                                ]
                                        },
                                        {
                                                "class_name": "TG-HDL",
                                                "contains_some":
                                                [
                                                        "TG",
                                                        "HDL"
                                                ]
                                        },
                                        {
                                                "class_name": "BMI",
                                                "contains_some":
                                                [
                                                        "BMI"
                                                ]
                                        }
                                ]
                        }
                }
        }
}
```

Using these rules and criteria, it is possible to define a wide variety of schema for placing
loci into custom, meaningfully different categories.

### `make_heatmaps`

The tool you've been waiting for! Visualize your colocalization results across dimensions in a heatmap using the 
locus numbers and categories you've defined in the previous steps.

Typically, each row of the heatmap will contain results for a single gene, with genes grouped logically
by locus number. Each column will correspond to a typical [tissue x trait] combination, with all
possible combinations represented in their own row. The colors of heatmap entries will then indicate
whether a particular combination of gene, tissue, and trait were tested, and whether they colocalized
for eQTLs, sQTLs, both, or neither. Currently the heatmap is tool is coded to recognize "eqtl" and "sqtl"
as valid QTL types for colocalization tests, although it could easily be modified to include other QTL types
and combinations as well. Similarly, although this format naturally lends itself to depicting the intersection
of tissues and traits on the x-axis, it would also be possible, for example, to substitute other natural QTL
contexts such as ancestry of the QTL individual in place of the QTL tissue field.

The required top-level parameters for generating a heatmap are as follows:
* `gwas_column`: The header entry for the column that indicates which trait was tested. If custom text is desired
	for plotting axis labels, this column can be defined using a mapping in the `mutate_columns` step.
* `tissue_column`: The header entry for the column that indicates which tissue was tested for QTLs.
* `type_column`: The header entry for the column that indicates which QTL type (eQTL or sQTL) was tested.

The following top-level parameters may be optionally specified to customize the output format.
* `tissue_order`: A list showing the order in which tissues should be displayed in the heatmap. These values 
	should correspond to the values contained in `tissue_column`.
* `gwas_order`: A list showing the order in which GWAS traits should be displayed in the heatmap. These values
        should correspond to the values contained in `gwas_column`.
* `cluster`: If `"True"`, perform hierarchical clustering to order the rows of the results heatmap rather than
	grouping them by locus. NOTE: This may not currently work as expected; further troubleshooting is needed.
* `put_scores_in_cells`: If `"True"`, put a numerical value in each cell of the heatmap indicating the colocalization
	score as a percentage. Works best for scores that are naturally interpretable as probabilities of colocalization,
	between 0 and 1.
* `rows_per_page`: If a heatmap contains more than this many rows, the additional rows will be extended on to another page.
* `y_axis_collapse`: If set to `"genes"`, then all genes at each locus are collapsed into a single row, with the new row
	entries representing the results at the locus level rather than the gene level. Default setting is "none".
* `x_axis_collapse`: If set to `"tissues"`, all columns (tissues) for each _GWAS trait_ will be collapsed into a single column,
	with the new column representing the results for the GWAS trait across all tissues. If set to `"gwas"`, all columns (GWAS traits)
	for each _tissue_ will be collapsed into a single column, with the new column representing the results for the tissue
	across all GWAS traits. If set to `"tissues-gwas"`, all columns across all tissues _and_ all GWAS traits will be collapsed
	into just one column, with the new column representing the gene's or locus's colocalization status across every tested
	combination of traits and tissues. Default setting is "none".
* `locus_selection_list`: Optionally, a list of locus numbers. The tool will output results only for locus numbers contained in this list.
* `gene_selection_list`: Optionally, a list of Ensembl gene IDs. The tool will output results only for genes contained in this list.

Finally, an additional parameter `file_strata` is required, containing a list of objects / dictionaries.
For each object, a complete set of output heatmaps will be generated, but the specifications with further
detail how each plot is to be arranged.

Each object in `file_strata` must contain the following parameters:
* `out_dir`: The output directory for the heatmaps, either as a path relative to the working directory
	or as an absolute path.
* `concise`: If `"True"`, remove all genes from the output that were tested but had no colocalizations.
	If `"False"`, all candidate genes overlapping GWAS loci will be included, even if they had no 
	colocalization under any tested condition.

Each object in `file_strata` can also optionally contain the parameters:
* `split_factors`: A list of column headers on which to sort the output heatmaps. For example,
	in the second of the three objects in the `file_strata` for the example config file below,
	loci in the results will be sorted by tissue-specificity, and an individual heatmap will be output
	containing each category of tissue-specific loci. Generally, this field will contain one of the
	columns created by the `classify_results` tool.
* `constrain_split`: Optionally, a list of allowable values for the column defined with `split_factors`.
	All other values in this column will be removed.
* `gwas_blacklist`: A list of values for `gwas_column` that should be excluded from the output heatmap.

#### Example

This example from `config/ir.config` creates three different sets of heatmaps: one arranged by the
number of candidate and colocalized genes, one arranged by tissue specificity, and one arranged by the
priority level of the highest-priority traits that colocalized for a locus.

```
{
        "make_heatmaps":
        {
                "gwas_column": "gwas_short_trait",
                "tissue_column": "qtl_tissue",
                "type_column": "qtl_type",
                "tissue_order": ["Subcutaneous-Adipose", "Visceral-Adipose",  "Skeletal-Muscle", "Liver", "Pancreas"],
                "gwas_order": ["MI", "ISI", "FI", "FG", "T2D","WHR","HDL","TG","BMI"],
                "cluster": "False",
                "put_scores_in_cells": "True",
                "y_axis_collapse": "none",
                "x_axis_collapse": "none",
                "file_strata":
                [
                        {
                                "out_dir": "concise/step1/",
                                "split_factors": ["step1_coloc_sorting"],
                                "concise": "True",
                                "gwas_blacklist": ["ISI"]
                        },
                        {
                                "out_dir": "concise/step2/",
                                "split_factors": ["step2_tissue_sorting"],
                                "concise": "True",
                                "gwas_blacklist": ["ISI"]
                        },
                        {
                                "out_dir": "concise/step3/",
                                "split_factors": ["step3_gwas_sorting"],
                                "concise": "True",
                                "gwas_blacklist": ["ISI"]
                        }
                ]
	}
}
```
