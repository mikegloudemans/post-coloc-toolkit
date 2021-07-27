import sys
import os
import numpy as np
import pprint
pp = pprint.PrettyPrinter(indent=4)

####################################################################
# post-colocalization Wizard!
####################################################################

# An interactive wizard for getting the toolkit up and running ASAP!

def main():

	#############################################
	# Part 0: Get input file
	#############################################

	message = '''\nWelcome to the post-colocalization toolkit!\n\nTo run these tools, you'll need a single tab-delimited file where each row represents the results of a colocalization test.\n\nIf you have such a file, enter the absolute or relative path here. You can type 'exit' at any point to leave the wizard and then come back later where you left off.'''

	# Set up config file
	config = {}

	valid = False
	inp = screen_input(message, config)
	while not valid:
		
		if inp.lower() == "exit":
			sys.exit()
		if not os.path.isfile(inp):
			print("\n\nThat file doesn't seem to exist. Enter a valid file name, or 'exit' to leave the wizard and double-check.\n\n")
		else:
			coloc_file = inp
			valid = True
			num_results = 0

			with open(coloc_file) as f:
				header = f.readline().strip().split("\t")
				for line in f:
					num_results += 1

			print("\n\nGreat, we found a file with the following columns:")
			print(header)
			confirm = screen_input(f"\nAnd it contains results from {num_results} colocalization tests. Does that sound right? (yes/no)\n\n", config)
			if not confirm.lower() == "yes":
				valid=False
				inp = screen_input("\n\nOK then, let's try again. Please enter the file path without quotation marks or extra symbols, or type 'exit' to leave the wizard.\n\n", config)

	# TODO: Do validation for all required columns here, maybe just in a loop actually
	
	required_columns = ["ref_snp", "qtl_file", "gwas_trait", "feature", "neg_log_gwas_pval", "neg_log_qtl_pval", "score", "ensembl"]
	for column in required_columns:

		if column not in header:
			print(f"The input file needs to have the following columns: {required_columns}. I didn't detect a '{column}' column. Please reformat the file accordingly and then re-enter the wizard.")
			sys.exit()
	
	# Now get the key info from the file
	scores = []
	with open(inp) as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			scores.append(data[header.index("score")])

	#############################################
	# Part 1: Post-hoc filter
	#############################################

	config["post_hoc_filter"] = {}

	message = print("Ready to configure step 1: post-hoc-filter...\n")

	print("First, let's determine a cutoff score for colocalization.\n\nFor reference, the quantiles of scores in your input data are the following:\n")
	quantiles = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
	score_quants = np.quantile(scores, quantiles)
	for i in range(len(score_quants)):
		print(f"{quantiles[i]*100}%:\t{score_quants[i]}")

	done = False
	while done != True:
		inp = screen_input("Please enter a cutoff score:\n", config)
		try:
			threshold = float(inp)
			done = True
		except:
			print("Whoops, that doesn't look right. The cutoff score should be a numerical value.\n")
			continue

	config["post_hoc_filter"]["colocalization_threshold"] = threshold

	yesno = get_yes_no("Excellent.\n\nI found the following GWAS traits in your file:\n\n{all_gwas_traits}\n\nOptionally, you can either exclude specific traits, or choose a specific limited set of traits to include.\n\nWould you like to exclude any specific traits? (yes/no)")
	if yesno:
		exclude_traits = set([])
		
		while True:
			inp = screen_input("Enter a trait to exclude, or type 'done' if you don't have any others to add. Or 'restart' to start from an empty list again.\n\n")

			if inp.lower() == "done":
				break

			if inp.lower() == "restart":
				continue

			found_trait = False
			for agt in all_gwas_traits:
				if agt.lower() == inp.lower()
					exclude_traits.add(agt)
					found_trait = True
					break
			
			if not found_trait:
				print(f"\n\nHmm, I didn't find that trait anywhere in the results.\n\nCurrently excluded traits:{list(exclude_traits)}\n")

		config["post_hoc_filter"]["removed_gwas"] = list(exclude_traits)

	yesno = get_yes_no("\n\nNow, do you want to specify a narrow list of traits to include, and remove all others? (yes/no)")

        if yesno:
                include_traits = set([])

                while True:
                        inp = screen_input("Enter a trait to include, or type 'done' if you don't have any others to add. Or 'restart' to start from an empty list again.\n\n")

                        if inp.lower() == "done":
                                break

                        if inp.lower() == "restart":
                                continue

                        found_trait = False
                        for agt in all_gwas_traits:
                                if agt.lower() == inp.lower()
                                        include_traits.add(agt)
                                        found_trait = True
                                        break

                        if not found_trait:
                                print(f"\n\nHmm, I didn't find that trait anywhere in the results.\n\nCurrently included traits:{list(include_traits)}\n")

                config["post_hoc_filter"]["kept_gwas"] = list(include_traits)

	yesno = get_yes_no("\n\nNext, I found the following QTL files in your input:\n\n{all_qtl_files}\n\nOptionally, you can either exclude specific files, or choose a specific limited set to include.\n\nWould you like to exclude any specific files? (yes/no)")

        if yesno:
                exclude_qtls = set([])

                while True:
                        inp = screen_input("Enter a QTL file to exclude, or type 'done' if you don't have any others to add. Or 'restart' to start from an empty list again.\n\n")

                        if inp.lower() == "done":
                                break

                        if inp.lower() == "restart":
                                continue

                        found_trait = False
                        for aqf in all_qtl_files:
                                if aqf.lower() == inp.lower()
                                        exclude_qtls.add(aqf)
                                        found_trait = True
                                        break

                        if not found_trait:
                                print(f"\n\nHmm, I didn't find that QTL file anywhere in the results.\n\nCurrently excluded:{list(exclude_qtls)}\n")

                config["post_hoc_filter"]["removed_qtl"] = list(exclude_qtls)

	yesno = get_yes_no("\n\nDo you want to specify a narrow list of QTL files to include, and remove all others? (yes/no)")

        if yesno:
                include_qtls = set([])

                while True:
                        inp = screen_input("Enter a QTL file to include, or type 'done' if you don't have any others to add. Or 'restart' to start from an empty list again.\n\n")

                        if inp.lower() == "done":
                                break

                        if inp.lower() == "restart":
                                continue

                        found_trait = False
                        for aqf in all_qtl_files:
                                if aqf.lower() == inp.lower()
                                        include_qtls.add(aqf)
                                        found_trait = True
                                        break

                        if not found_trait:
                                print(f"\n\nHmm, I didn't find that QTL file anywhere in the results.\n\nCurrently included:{list(include_qtls)}\n")

                config["post_hoc_filter"]["kept_qtl"] = list(include_qtls)

	print("\n\nYou can also apply filters to the p-values for GWAS and QTLs. If you want to set different p-value thresholds for different GWAS traits or for different QTL files, please see the README documentation and modify the resulting config file accordingly. Here, we'll assume one p-value threshold for all the GWAS traits and one for all the QTL files.\n")
	yesno = get_yes_no("Would you like to filter results by GWAS p-value? (yes/no)")
	
	if yesno:

		print("\n\nEnter the GWAS p-value threshold (unadjusted, not on a log scale). You can use scientific notation if you'd like, e.g. '5e-8'. Or type 'skip' if you've decided not to include a GWAS p-value threshold.\n\n")

		no_threshold = False

		while True:

			inp = screen_input("Enter a value:\n\n")

			if inp.lower() == 'skip':
				no_threshold = True
				break

			try:
				gwas_pval_threshold = float(inp)
				if gwas_pval_threshold < 0 or gwas_pval_threshold > 1: 
					print("\n\nThat doesn't fall within the p-value range of 0 to 1 a numerical value. Try again.\n")
					continue
				break
			except:
				print("\n\nThat doesn't look like a numerical value. Try again.\n")
				continue

		if not no_threshold:
			config["post_hoc_filter"]["gwas_pval_threshold"]= {"standard": gwas_pval_threshold}

			print("Saved GWAS p-value threshold.\n\n")


	yesno = get_yes_no("Would you like to filter results by QTL p-value? (yes/no)")

	if yesno:

		print("\n\nEnter the QTL p-value threshold (unadjusted, not on a log scale). You can use scientific notation if you'd like, e.g. '5e-8'. Or type 'skip' if you've decided not to include a QTL p-value threshold.\n\n")

		no_threshold = False

		while True:

			inp = screen_input("Enter a value:\n\n")

			if inp.lower() == 'skip':
				no_threshold = True
				break

			try:
				qtl_pval_threshold = float(inp)
				if qtl_pval_threshold < 0 or qtl_pval_threshold > 1: 
					print("\n\nThat doesn't fall within the p-value range of 0 to 1 a numerical value. Try again.\n")
					continue
				break
			except:
				print("\n\nThat doesn't look like a numerical value. Try again.\n")
				continue

		if not no_threshold:
			config["post_hoc_filter"]["qtl_pval_threshold"]= {"standard": qtl_pval_threshold}

			print("Saved QTL p-value threshold.\n\n")

	yesno = get_yes_no("Finally, would you like to set a minimum number of SNPs in a test for it to be counted?")

	if yesno:
	
		no_threshold = False

		while True:

			inp = screen_input("Enter the minimum number of allowed SNPs. (If you're not sure, 20 is a good starting point.) Or type 'skip' if you've decided not to set a minimum threshold.\n\n")

			if inp.lower() == 'skip':
				no_threshold = True
				break

			try:
				snps = int(inp)
				if snp_threshold < 1: 
					print("\n\nNice try...you need a positive integer number of SNPs.\n")
					continue
				break
			except:
				print("\n\nNice try...you need a positive integer number of SNPs.\n")
				continue

		if not no_threshold:
			config["post_hoc_filter"]["snp_count_min_threshold"]= snps

			print("Saved minimum SNP count threshold.\n\n")

	config["post_hoc_filter"]["completed_wizard"] = "True"	

	print("Awesome, we're through configuring the post_hoc_filter tool.\n"
	
	#############################################
	# Part 2: Mutate columns
	#############################################

	config["mutate_columns"] = {"mutations": []}

	print("Ready to configure step 2: mutate_columns...\n")
	
	print("This step will create new columns from the existing column names (while preserving the original column names). We'll walk through a few examples where you might want to do this, and then you can add an others you might want.")

	yesno = get_yes_no(f"The GWAS traits in your file are {all_gwas_traits}. Would you like to map these to a new column? You might want to do this if you want to shorten the names for display, if multiple studies correspond to the same trait, or if the traits naturally fall in to multiple categories. (You can do multiple remappings if needed.)\n\nAnswer (yes/no):")
	while yesno:

		out_column = screen_input("\n\nWhat do you want the name of your new column to be? (e.g. 'short_gwas_name')\n\n")
		
		mapping = {"in": "gwas_trait", "out": out_column, "map": {}}

		print("For each field in the column, enter the remapped field you'd like to place in the new column...")

		for field in all_gwas_traits:
			out_field = screen_input(f"{field}: ")
			mapping["map"][field] = out_field

		print("Great, this is the mapping you specified...")
		pp.pprint(mapping)
		nested_yesno = get_yes_no("Does this look right? (yes/no)")

		if not nested_yesno:
			print("OK, you can try specifying it again. We'll drop this mapping from the config file.")
		else:
			config["mutate_columns"]["mutations"].append(mapping)
			print(f"Added column remapping 'gwas_trait':'{out_column}' to the list.")

		yesno = get_yes_no("Do you want to add any other remappings for the GWAS traits? (yes/no)")
	
	yesno = get_yes_no(f"The QTL files in your file are {all_qtl_files}. Would you like to map these to a new column? You might want to do this if you want to shorten the names for display, if you have multiple files for the same tissue, or if you have QTL types (e.g. eQTL and sQTL) that you'd like to specify. (You can do multiple remappings if needed.)\n\nAnswer (yes/no):")
	while yesno:

		out_column = screen_input("\n\nWhat do you want the name of your new column to be? (e.g. 'qtl_tissue')\n\n")
		
		mapping = {"in": "qtl_file", "out": out_column, "map": {}}

		print("For each field in the column, enter the remapped field you'd like to place in the new column...")

		for field in all_qtl_files:
			out_field = screen_input(f"{field}: ")
			mapping["map"][field] = out_field

		print("Great, this is the mapping you specified...")
		pp.pprint(mapping)
		nested_yesno = get_yes_no("Does this look right? (yes/no)")

		if not nested_yesno:
			print("OK, you can try specifying it again. We'll drop this mapping from the config file.")
		else:
			config["mutate_columns"]["mutations"].append(mapping)
			print(f"Added column remapping 'qtl_file':'{out_column}' to the list.")

		yesno = get_yes_no("Do you want to add any other remappings for the QTL files? (yes/no)")
	
	yesno = get_yes_no(f"The rest of the columns in your file are {header}. Do you want to remap any of these other columns? (For most applications, the answer is probably 'no'). (yes/no)")
	while yesno:

		in_column = screen_input("\n\nWhich input column do you want to remap?\n\n")

                found_trait = False
                for h in header:
                        if h.lower() == in_column.lower()
                                found_trait = True
                                break

                if not found_trait:
                        yesno = get_yes_no(f"\n\nHmm, that doesn't seem to be one of the columns in the header. Do you still want to remap one of the columns? (yes/no)\n")
			continue

		out_column = screen_input("\n\nWhat do you want the name of your new column to be? (e.g. 'qtl_tissue')\n\n")
		
		mapping = {"in": in_column, "out": out_column, "map": {}}

		print("For each field in the column, enter the remapped field you'd like to place in the new column...")

		for field in # TODO: we need to get the fields now...:
			out_field = screen_input(f"{field}: ")
			mapping["map"][field] = out_field

		print("Great, this is the mapping you specified...")
		pp.pprint(mapping)
		nested_yesno = get_yes_no("Does this look right? (yes/no)")

		if not nested_yesno:
			print("OK, you can try specifying it again. We'll drop this mapping from the config file.")
		else:
			config["mutate_columns"]["mutations"].append(mapping)
			print(f"Added column remapping 'qtl_file':'{out_column}' to the list.")

		yesno = get_yes_no("Do you want to add any other remappings for the other columns? (yes/no)")
	
	config["mutate_columns"]["completed_wizard"] = "True"	

	print("Awesome, we're through configuring the mutate_columns tool.\n"

	#############################################
	# Part 3: Add HGNC naems
	#############################################

def get_yes_no(message, part_config):

	inp = screen_input(message, part_config)

	while inp.lower() not in ["yes", "no"]:
	
		message = "Please enter 'yes' or 'no'.\n\n"

		inp = screen_input(message, part_config)

	return(inp == "yes")

def screen_input(message, part_config):

	inp = input(message).strip()

	if inp.lower() == "exit":
		if len(part_config.keys()) != 0:
			# TODO: Write temporary config files
			print("OK, exiting wizard and saving your progress.")
			pass
		sys.exit()

	return(inp)

main()

