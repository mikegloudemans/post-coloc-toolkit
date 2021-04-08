import sys

####################################################################
# post-colocalization Wizard!
####################################################################

# An interactive wizard for getting the toolkit up and running ASAP!

def main()
    message = '''
    Welcome to the post-colocalization toolkit!\n\nTo run this kit, you'll
    need a single file where each row represents the results of a colocalization test.\n\n
    If you have such a file, enter the absolute or relative path here, 
    or type "exit" to leave the wizard and prepare
    your file.'''

    valid = False
    while not valid:
        inp = input(message)
        valid, message = validate_results_file(inp)

def validate_results_file(inp)
    if inp == "exit":
        sys.exit()
    if not os.path.isfile(inp):
        return False, "That file doesn't seem to exist. Enter a valid file name,
                or 'exit' to leave the wizard and double-check."
    return valid, None

main()

Wizard:
Thanks for trying out the post-colocalization-toolkit! We have just a few questions to get you started. Also, at any point, you can terminate this wizard part-way through by typing “exit” and pressing ENTER, and then when you start back up again, we’ll give the option of picking up where you left off.
First, you’ll need one or more tab-separated files containing your colocalization results. They also need to be named with the suffix “TODO”. Do you have one or more such files created?
NO -> all right, please create or rename those files and then come back when you’re ready. (END)
Please type the location of the folder in which these files are located and then press ENTER. It’s OK if your files are a few folders deep, but they need to be somewhere nested in this folder.:
None found: We didn’t find any colocalization files in this directory with the suffix “TODO”. Could you try re-entering the folder name?
	Some found: Great, we found (N) files (display the first 10 or so). Does this seem right?
NO: OK, maybe that wasn’t the right folder. Please type the name of the folder where the files are located. (back to that step)
Great, we’ve got our list of colocalization files, and we’ll combine these together shortly. Can you verify that these files are all in the same format (i.e. have the same number and order of columns, with the same column headers?)
	Yep, got it
	Wait, I have to fix this
	Awesome! We’ve combined these (X) files into 1 individual file with colocalization results for (Y) loci. Would you like us to save this intermediate file? (get location too)
		Or “Hmm…we tried combining those files together, but they didn’t fit…” (display message telling why)
	The next step is applying any post-hoc-filters you might be interested in including.
		We see that the GWAS files you’ve included are (list them). Would you like to remove some of these? 
			OR “we can’t find a ‘gwas_trait’ column in your table. This will be a problem if you want to filter out certain GWAS traits, so we’re going to skip this step. Is this OK?”
		(Same for QTL types)
		What about p-values? Your current range of GWAS p-values is from [X] to [Y]. Should we make this filter more stringent?
			Should this filter apply to just one individual study, or to all GWAS?
		OK, great. Your current GWAS filters are: (list them). Would you like to create any additional filters?
		Same for process for eQTL
		All right, that covers the filters. Just to confirm (list filters…) and restart if these filters aren’t right.
		Excellent. We’ve applied these filters and you now have X of Y colocalization tests remaining.

		We can run a few quick QC tests to make sure all the tests you wanted to run were properly executed. We highly recommend completing this step. Would you like to run the interactive QC step?
			[make interactive QC step]
		If one of the columns in your table contains Ensembl gene IDs, we can convert these to more readable gene names (e.g. ENSG0000004061 becomes “ABC1”). That column would need to have the header “qtl_feature”. 
		We’ll need to download a mapping file from Ensembl IDs to HGNC gene ids. (approximate size: X MB) Is it OK to do this and to place the file at `location`?
		We can also add rsids to your table if they’re not already present. (e.g. 1_191292 -> rs19294). This will be helpful for looking up your SNPs on the UCSC Genome Browser and in other functional databases. You will need a column “ref_snp” (describe format). Would you like to obtain rsids?
			
		We’ll need to download a dbSNP mapping file from chromosomal coordinates to rsids. (approximate size: X MB) Is it OK to do this and to place the file at `location`?


		Great! There’s one required annotation step, which is clustering nearby SNPs into the same locus. The default setting is to consider any lead SNPs less than 1 Mb apart as members of the same locus. This usually works well, but if you have a LOT of GWAS or eQTL studies included, this window size may be too big. Would you like to specify a custom locus size?
		OK, how far apart do two lead GWAS SNPs have to be for us to consider them separate loci? (A word of warning is that if you choose too small a number, you might end up calling two or more colocalization tests different loci, when they’re actually just tagging the same SNP.)

		OK, excellent. We’ve determined that your set of colocalization tests includes a total of N independent loci between the M total reference SNPs.
		(include quantile tool for exploring this…)
		Now for a more interesting step. We can categorize all the loci you’ve created in several different ways. We’ll first suggest a few common “rules” for categorization, and then you can create some custom rules later if you’d like.
		One more thing: what threshold would you like to use as the cutoff between “colocalization” and “no colocalization”? We can’t recommend a specific threshold because this largely depends on which method you’re using, but we can show you the percent of tests passing at each threshold in your set: 

		First of all, we might want to categorize loci based on whether they had any likely colocalizations, and if so, how many they had. A sample rule set is the following:
		Class 1: Single candidate gene w/ GWAS and eQTL, no colocalization
		Class 2: Single candidate gene w/ GWAS and eQTL, colocalization
		Class 3: Multiple candidates, no colocalization
		Class 4: Multiple candidates, multiple colocalizations
		Class 5: Multiple candidates, one colocalization
		Would you like us to include this classification schema as a default rule?

		Another rule we might use for classification is a QTL/GWAS type-based categorization. We see that you have X types of QTLs remaining. Thus, we could categorize our loci into those colocalizing in exactly one QTL type, those colocalizing in none, and those colocalizing in more than one QTL type. Should we add such a rule?
			If yes…
		All right, we added this rule. One quick word of warning: just because you see colocalization in a particular QTL type doesn’t guarantee that the effect is specific to that QTL type. It might just be that you didn’t test that locus in all QTL types, or that the relevant QTL type wasn’t tested, among other possibilities.

		One other type of rule we could use is a tier-based rule for prioritizing certain GWAS. We see that your results include N different types of GWAS. A tier-based rule would prioritize all SNPs colocalized in a certain GWAS…then give lower priority to all those prioritized in another GWAS but not in the first…and so on. Would you like to create a tier-based rule for these GWAS?
			If yes:
			OK, which GWAS would you consider the most important? (you can enter the number of the GWAS of interest.)
				Great, that GWAS is now your top tier.
				Would you like to add another tier?
					(repeat as needed)
				Excellent, we’ve now added a tier-based rule.

				That’s enough rules for now, but the possibilities for categorization go well beyond what you’ve seen here. See the documentation on Github accompanying the “classify_results” module for more information on how to finely customize this step.
				(what if they made no rules? We’ll just include an example then I guess, the colocalization-based one, but not put it in the config file for them.) 
				We ran the categorization with your desired rules, and obtained the following results:
				[print results]
				We’ve also added the categorizations for each coloc test to the table of colocalization results.

				That wraps up the direct characterization of the colocalization loci. However, there are a few other follow-up tests we can explore.

				Would you like to get a summary of all the filters performed at various steps within your colocalization run? This analysis can also optionally produce a graphic illustrating how many SNPs, genes, and loci were removed from consideration at various steps of the process.
					Great. We’ll do all the default filters.
				Since you ran the classify_results step, we can summarize an additional filter if you’re interested in filtering your results based on this column. Is this something you’d be interested in doing?
				(then just give options and whatever).
				Finally, would you like to make a graphic showing the different steps of filtering and how many SNPs, genes, and loci passed at each?

				One final analysis we can do is to get LD buddies of colocalized variants and even predict the potential consequences for each of these variants using VEP. Are you interested in either of these steps?
				[Get LD buddies]
				[Get VEP consequences]

				(By this point we’ll have run the whole pipeline, but we’ll also show a config file to clarify how this was all done and so it can be rerun as needed). 
				(If “break” is typed at any point, save partial config if needed and load from there so progress is not lost. Or give the option to start over, when we restart this wizard. Can also save file to a custom location at the end.)

