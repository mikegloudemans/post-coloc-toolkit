import sys
import os
import numpy as np
import pprint
import json
pp = pprint.PrettyPrinter(indent=4)

def save_results(part_config):
	with open(config_file, "w") as w:
		json.dump(part_config, w, indent=4)

def get_yes_no(message, part_config):

	inp = screen_input(message, part_config)

	while inp.lower() not in ["yes", "no"]:
	
		message = "\nPlease enter 'yes' or 'no'.\n\n"

		inp = screen_input(message, part_config)

	return(inp == "yes")

def screen_input(message, part_config):

	inp = input(message+">>  ").strip()

	if inp.lower() == "exit":
		if len(part_config.keys()) != 0:

			if "post_hoc_filter" not in config or "completed_wizard" not in config["post_hoc_filter"] or config["post_hoc_filter"]["completed_wizard"] != "True":
				print("OK, exiting wizard.")
				sys.exit()

			with open(config_file, "w") as w:
				json.dump(part_config, w, indent=4)

			# TODO: Write temporary config files
			print("OK, exiting wizard and saving your progress.")
			pass
		sys.exit()

	return(inp)

default_num_colocs_rule = {
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
			}


####################################################################
# post-colocalization Wizard!
####################################################################

# An interactive wizard for getting the toolkit up and running ASAP!

print("\nWelcome to the post-colocalization toolkit!\n")

config = {}
config_file = screen_input("Please enter the path where you'd like to create a config file. (e.g. for 'config/my-tests.config', just enter 'my-tests.config')\n\n", {})

config_file = f"config/{config_file}"

if os.path.isfile(config_file):
	yesno = get_yes_no("\nThat config file already exists. Would you like to load it and resume creating it? Otherwise, this file will be overwritten. (yes/no)\n\n", {})

	if yesno:
		with open(config_file) as f:
			config = json.load(f)

def main():

	#############################################
	# Part 0: Get input file
	#############################################

	message = '''\nTo run these tools, you'll need a single tab-delimited file where each row represents the results of a colocalization test.\n\nIf you have such a file, place it in 'data/coloc_results/' and then enter the actual filename here. You can type 'exit' at any point to leave the wizard and then come back later where you left off.\n\n'''

	valid = False
	inp = screen_input(message, config)

	while True:
		
		if inp.lower() == "exit":
			sys.exit()
		if not os.path.isfile(f"data/coloc_results/{inp}"):
			inp = screen_input("\n\nThat file doesn't seem to exist. Enter a valid file name, or 'exit' to leave the wizard and double-check.\n\n", config)
			continue
		else:
			coloc_file = f"data/coloc_results/{inp}"
			num_results = 0

			with open(coloc_file) as f:
				header = f.readline().strip().split("\t")
				for line in f:
					num_results += 1

			print("\nGreat, we found a file with the following columns:")
			print(header)
			confirm = screen_input(f"\nAnd it contains results from {num_results} colocalization tests. Does that sound right? (yes/no)\n\n", config)
			if not confirm.lower() == "yes":
				inp = screen_input("\n\nOK then, let's try again. Please enter the file path without quotation marks or extra symbols, or type 'exit' to leave the wizard.\n\n", config)
				continue

			if config_file.replace("config/", "").replace(".config", "") != coloc_file.replace("data/coloc_results/", "").replace("_colocalization_results.txt", ""):
				print(f"\nALERT: You've currently specified {config_file} as your config file and {coloc_file} as your input file. For the pipeline to run, you'll eventually need them to have the same prefix (e.g. 'ir_colocalization_results.txt' and 'ir.config'). You can rename one of the files later to fix this.\n")

			break

	# TODO: Do validation for all required columns here, maybe just in a loop actually
	
	required_columns = ["ref_snp", "qtl_file", "gwas_trait", "feature", "neg_log_gwas_pval", "neg_log_qtl_pval", "score", "ensembl"]
	for column in required_columns:

		if column not in header:
			print(f"The input file needs to have the following columns: {required_columns}. I didn't detect a '{column}' column. Please reformat the file accordingly and then re-enter the wizard.\n\n")
			sys.exit()
	
	# Now get the key info from the file
	scores = []
	all_gwas_traits = set([])
	all_qtl_files = set([])
	with open(coloc_file) as f:
		f.readline()
		for line in f:
			data = line.strip().split()
			scores.append(float(data[header.index("score")]))
			all_gwas_traits.add(data[header.index("gwas_trait")])
			all_qtl_files.add(data[header.index("qtl_file")])

	if "post_hoc_filter" not in config or "completed_wizard" not in config["post_hoc_filter"] or config["post_hoc_filter"]["completed_wizard"] != "True":

		#############################################
		# Part 1: Post-hoc filter
		#############################################

		config["post_hoc_filter"] = {}

		message = print("\nReady to configure step 1: post-hoc-filter...\n")

		print("First, let's determine a cutoff score for colocalization.\n\nFor reference, the quantiles of scores in your input data are the following:\n")
		quantiles = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
		score_quants = np.quantile(scores, quantiles)
		for i in range(len(score_quants)):
			print(f"{quantiles[i]*100}%:\t{score_quants[i]}")

		done = False
		while done != True:
			inp = screen_input("\nPlease enter a cutoff score:\n\n", config)
			try:
				threshold = float(inp)
				done = True
			except:
				print("\nWhoops, that doesn't look right. The cutoff score should be a numerical value.\n")
				continue

		config["post_hoc_filter"]["colocalization_threshold"] = threshold

		yesno = get_yes_no(f"\nExcellent.\n\nI found the following GWAS traits in your file:\n\n{all_gwas_traits}\n\nOptionally, you can either exclude specific traits, or choose a specific limited set of traits to include.\n\nWould you like to exclude any specific traits? (yes/no)\n\n", config)
		if yesno:
			exclude_traits = set([])
			
			while True:

				print(f"\nCurrently excluded traits: {exclude_traits}")

				inp = screen_input("\nEnter a trait to exclude, or type 'done' if you don't have any others to add. Or 'restart' to start from an empty list again.\n\n", config)

				if inp.lower() == "done":
					break

				if inp.lower() == "restart":
					exclude_traits = set([])
					continue

				if inp not in all_gwas_traits:
					print(f"\nHmm, I didn't find that trait anywhere in the results.")
					continue

				exclude_traits.add(inp)

			config["post_hoc_filter"]["removed_gwas"] = list(exclude_traits)

		yesno = get_yes_no("\nNow, do you want to specify a narrow list of traits to include, and remove all others? (yes/no)\n\n", config)

		if yesno:
			include_traits = set([])

			while True:
				print(f"\nCurrently included traits: {include_traits}")

				inp = screen_input("\nEnter a trait to include, or type 'done' if you don't have any others to add. Or 'restart' to start from an empty list again.\n\n", config)

				if inp.lower() == "done":
					break

				if inp.lower() == "restart":
					include_traits = set([])
					continue

				if inp not in all_gwas_traits:
					print(f"\nHmm, I didn't find that trait anywhere in the results.")
					continue
				
				include_traits.add(inp)

			config["post_hoc_filter"]["kept_gwas"] = list(include_traits)

		yesno = get_yes_no(f"\n\nNext, I found the following QTL files in your input:\n\n{all_qtl_files}\n\nOptionally, you can either exclude specific files, or choose a specific limited set to include.\n\nWould you like to exclude any specific files? (yes/no)\n\n", config)

		if yesno:
			exclude_qtls = set([])

			while True:
				print(f"\nCurrently excluded files: {exclude_qtls}")

				inp = screen_input("Enter a QTL file to exclude, or type 'done' if you don't have any others to add. Or 'restart' to start from an empty list again.\n\n", config)

				if inp.lower() == "done":
					break

				if inp.lower() == "restart":
					continue

				found_trait = False
				for aqf in all_qtl_files:
					if aqf.lower() == inp.lower():
						exclude_qtls.add(aqf)
						found_trait = True
						break

				if not found_trait:
					print(f"\n\nHmm, I didn't find that QTL file anywhere in the results.\n\nCurrently excluded:{list(exclude_qtls)}\n")

			config["post_hoc_filter"]["removed_qtl"] = list(exclude_qtls)

		yesno = get_yes_no("\nDo you want to specify a narrow list of QTL files to include, and remove all others? (yes/no)\n\n", config)

		if yesno:
			include_qtls = set([])

			while True:
				print(f"\nCurrently included files: {include_qtls}")

				inp = screen_input("Enter a QTL file to include, or type 'done' if you don't have any others to add. Or 'restart' to start from an empty list again.\n\n", config)

				if inp.lower() == "done":
					break

				if inp.lower() == "restart":
					continue

				found_trait = False
				for aqf in all_qtl_files:
					if aqf.lower() == inp.lower():
						include_qtls.add(aqf)
						found_trait = True
						break

				if not found_trait:
					print(f"\n\nHmm, I didn't find that QTL file anywhere in the results.\n\nCurrently included:{list(include_qtls)}\n")

			config["post_hoc_filter"]["kept_qtl"] = list(include_qtls)

		print("\nYou can also apply filters to the p-values for GWAS and QTLs. If you want to set different p-value thresholds for different GWAS traits or for different QTL files, please see the README documentation and modify the resulting config file accordingly. Here, we'll assume one p-value threshold for all the GWAS traits and one for all the QTL files.\n")
		yesno = get_yes_no("Would you like to filter results by GWAS p-value? (yes/no)\n\n", config)
		
		if yesno:

			print("\nEnter the GWAS p-value threshold (unadjusted, not on a log scale). You can use scientific notation if you'd like, e.g. '5e-8'. Or type 'skip' if you've decided not to include a GWAS p-value threshold.\n")

			no_threshold = False

			while True:

				inp = screen_input("Enter a value:\n\n", config)

				if inp.lower() == 'skip':
					no_threshold = True
					break

				try:
					gwas_pval_threshold = float(inp)
					if gwas_pval_threshold < 0 or gwas_pval_threshold > 1: 
						print("\nThat doesn't fall within the p-value range of 0 to 1 a numerical value. Try again.\n")
						continue
					break
				except:
					print("\nThat doesn't look like a numerical value. Try again.\n")
					continue

			if not no_threshold:
				config["post_hoc_filter"]["gwas_pval_threshold"]= {"standard": gwas_pval_threshold}

				print("\nSaved GWAS p-value threshold.\n")


		yesno = get_yes_no("Would you like to filter results by QTL p-value? (yes/no)\n\n", config)

		if yesno:

			print("\n\nEnter the QTL p-value threshold (unadjusted, not on a log scale). You can use scientific notation if you'd like, e.g. '5e-8'. Or type 'skip' if you've decided not to include a QTL p-value threshold.\n")

			no_threshold = False

			while True:

				inp = screen_input("Enter a value:\n\n", config)

				if inp.lower() == 'skip':
					no_threshold = True
					break

				try:
					qtl_pval_threshold = float(inp)
					if qtl_pval_threshold < 0 or qtl_pval_threshold > 1: 
						print("\nThat doesn't fall within the p-value range of 0 to 1 a numerical value. Try again.\n")
						continue
					break
				except:
					print("\nThat doesn't look like a numerical value. Try again.\n")
					continue

			if not no_threshold:
				config["post_hoc_filter"]["qtl_pval_threshold"]= {"standard": qtl_pval_threshold}

				print("\nSaved QTL p-value threshold.\n")

		yesno = get_yes_no("Finally, would you like to set a minimum number of SNPs in a test for it to be counted?\n\n", config)

		if yesno:
		
			no_threshold = False

			while True:

				inp = screen_input("Enter the minimum number of allowed SNPs. (If you're not sure, 20 is a good starting point.) Or type 'skip' if you've decided not to set a minimum threshold.\n\n", config)

				if inp.lower() == 'skip':
					no_threshold = True
					break

				try:
					snps = int(inp)
					if snps < 1: 
						print("\nNice try...you need a positive integer number of SNPs.\n")
						continue
					break
				except:
					print("\nNice try...you need a positive integer number of SNPs.\n")
					continue

			if not no_threshold:
				config["post_hoc_filter"]["snp_count_min_threshold"]= snps

				print("Saved minimum SNP count threshold.\n\n")

		config["post_hoc_filter"]["completed_wizard"] = "True"	

		print("Awesome, we're through configuring the post_hoc_filter tool.\n")

		save_results(config)	
	
	else:	
		print("Step 1: post_hoc_filter has already been configured...")

	if "mutate_columns" not in config or "completed_wizard" not in config["mutate_columns"] or config["mutate_columns"]["completed_wizard"] != "True":
		#############################################
		# Part 2: Mutate columns
		#############################################

		config["mutate_columns"] = {"mutations": []}

		print("\nReady to configure step 2: mutate_columns...\n")
		
		print("This step will create new columns from the existing column names (while preserving the original column names). We'll walk through a few examples where you might want to do this, and then you can add an others you might want.\n")

		yesno = get_yes_no(f"The GWAS traits in your file are {all_gwas_traits}.\n\nWould you like to map these to a new column? You might want to do this if you want to shorten the names for display, if multiple studies correspond to the same trait, or if the traits naturally fall into multiple categories. (You can do multiple remappings if needed.)\n\nAnswer (yes/no):\n\n", config)
		while yesno:

			out_column = screen_input("\n\nWhat do you want the name of your new column to be? (e.g. 'short_gwas_name')\n\n", config)
			
			mapping = {"in": "gwas_trait", "out": out_column, "map": {}}

			print("For each field in the column, enter the remapped field you'd like to place in the new column. (or just press enter to keep the old value for that field).\n")

			for field in all_gwas_traits:
				out_field = screen_input(f"{field}: ", config)
				if out_field == "":
					mapping["map"][field] = field
				else:
					mapping["map"][field] = out_field			

			print("Great, this is the mapping you specified...")
			pp.pprint(mapping)
			print("\n")
			nested_yesno = get_yes_no("Does this look right? (yes/no)\n\n", config)

			if not nested_yesno:
				print("\nOK, you can try specifying it again. We'll drop this mapping from the config file.\n")
			else:
				config["mutate_columns"]["mutations"].append(mapping)
				print(f"\nAdded column remapping 'gwas_trait':'{out_column}' to the list.\n")

			yesno = get_yes_no("Do you want to add any other remappings for the GWAS traits? (yes/no)\n\n", config)
		
		yesno = get_yes_no(f"\nThe QTL files in your file are {all_qtl_files}.\n\nWould you like to map these to a new column? You might want to do this if you want to shorten the names for display, if you have multiple files for the same tissue, or if you have QTL types (e.g. eQTL and sQTL) that you'd like to specify. (You can do multiple remappings if needed.)\n\nAnswer (yes/no):\n\n", config)
		while yesno:

			out_column = screen_input("\n\nWhat do you want the name of your new column to be? (e.g. 'qtl_tissue')\n\n", config)
			
			mapping = {"in": "qtl_file", "out": out_column, "map": {}}

			print("For each field in the column, enter the remapped field you'd like to place in the new column. (or just press enter to keep the old value for that field).\n")

			for field in all_qtl_files:
				out_field = screen_input(f"{field}: ", config)
				if out_field == "":
					mapping["map"][field] = field
				else:
					mapping["map"][field] = out_field

			print("Great, this is the mapping you specified...")
			pp.pprint(mapping)
			nested_yesno = get_yes_no("Does this look right? (yes/no)\n\n", config)

			if not nested_yesno:
				print("OK, you can try specifying it again. We'll drop this mapping from the config file.")
			else:
				config["mutate_columns"]["mutations"].append(mapping)
				print(f"\nAdded column remapping 'qtl_file':'{out_column}' to the list.\n")

			yesno = get_yes_no("Do you want to add any other remappings for the QTL files? (yes/no)\n\n", config)
		
		yesno = get_yes_no(f"\nThe rest of the columns in your file are {header}. Do you want to remap any of these other columns? (For most applications, the answer is probably 'no'). (yes/no)\n\n", config)
		while yesno:

			in_column = screen_input("\n\nWhich input column do you want to remap?\n\n", config)

			if not in_column in header:
				yesno = get_yes_no(f"\n\nHmm, that doesn't seem to be one of the columns in the header. Do you still want to remap one of the columns? (yes/no)\n", config)
				continue

			out_column = screen_input("\n\nWhat do you want the name of your new column to be? (e.g. 'qtl_tissue')\n\n", config)
			
			mapping = {"in": in_column, "out": out_column, "map": {}}

			print("For each field in the column, enter the remapped field you'd like to place in the new column...or just press Enter to keep the same value for that field.\n")

			all_fields = set([])
			with open(coloc_file) as f:
				head = f.readline()
				h = head.index(in_column)
				for line in f:
					data = line.strip().split()
					all_fields.add(data[h])

			for field in list(all_fields):
				out_field = screen_input(f"{field}: ", config)
				mapping["map"][field] = out_field

			print("Great, this is the mapping you specified...")
			pp.pprint(mapping)
			nested_yesno = get_yes_no("Does this look right? (yes/no)", config)

			if not nested_yesno:
				print("OK, you can try specifying it again. We'll drop this mapping from the config file.")
			else:
				config["mutate_columns"]["mutations"].append(mapping)
				print(f"\nAdded column remapping 'qtl_file':'{out_column}' to the list.\n")

			yesno = get_yes_no("Do you want to add any other remappings for the other columns? (yes/no)\n\n", config)
		
		config["mutate_columns"]["completed_wizard"] = "True"	

		print("Awesome, we're through configuring the mutate_columns tool.\n")
		
		save_results(config)	
	
	else:	
		print("Step 2: mutate_columns has already been configured...")
	
	if "add_hgnc_names" not in config or "completed_wizard" not in config["add_hgnc_names"] or config["add_hgnc_names"]["completed_wizard"] != "True":

		#############################################
		# Part 3: Add HGNC names
		#############################################

		print("Now we'll configure step 3: add_hgnc_names...\n")

		config["add_hgnc_names"] = {}

		yesno = get_yes_no("A default file for mapping Ensembl gene IDs to HGNC gene ids is located at 'data/hgnc/ensembl_to_hgnc.txt'. You can specify a different mapping file if you want, though. Do you want to use a different mapping file? (yes/no)\n\n", config)

		header_yesno = "True"
		ensembl_col = 1
		hgnc_col = 3
		if yesno:
			while True:
				inp = screen_input("\nType the location of the file you'd like to use.\n\n", config)
		
				if not os.path.isfile(inp):
					print("That file doesn't seem to exist, please try again. Be sure either to use an absolute path or to specify the path relative to the 'wizard.py' directory.\n\n")
					continue

				print("Great. The head of your file is as follows:\n")
				with open(inp) as f:
					for i in range(10):
						print(f.readline().strip())

				while True:
					ensembl_col = screen_input("\nPlease enter the 1-based index of the column containing the Ensembl gene IDs.\n\n", config)
					try:
						ensembl_col = int(ensembl_col)
						break
					except:
						print("\nYou'll need an integer value for this parameter.")

				while True:
					hgnc_col = screen_input("\nPlease enter the 1-based index of the column containing the HGNC gene IDs.\n\n", config)
					try:
						hgnc_col = int(hgnc_col)
						break
					except:
						print("You'll need an integer value for this parameter.\n")

				header_yesno = get_yes_no("\nGreat. One last thing: does your file have a header row? (yes/no)\n\n", config)

				print("OK, so your mapping looks something like this...")

				map_head = {}
				with open(inp) as f:
					if header_yesno:
						f.readline()
					for i in range(10):
						data = f.readline().strip().split()
						ensembl = data[ensembl_col-1]
						hgnc = data[hgnc_col-1]
						map_head[ensembl] = hgnc

				pp.pprint(map_head)

				done = get_yes_no("Does that look right? (yes/no)\n\n", config)

				if not done:
					persist = get_yes_no("\nOops. We'll just drop this mapping then. Do you still want to specify an alternate gene mapping file? (yes/no)\n\n", config)
					if persist:
						continue

				break

		config["add_hgnc_names"]["header"] = header_yesno
		config["add_hgnc_names"]["ensembl_col_index"] = ensembl_col
		config["add_hgnc_names"]["hgnc_col_index"] = hgnc_col
		config["add_hgnc_names"]["ensembl_to_hgnc_map_file"] = inp

		config["add_hgnc_names"]["completed_wizard"] = "True"	

		print("\nAnd that's all we need for the add_hgnc_names tool.\n")
		
		save_results(config)	
	
	else:	
		print("Step 3: add_hgnc_names has already been configured...")

	if "assign_locus_numbers" not in config or "completed_wizard" not in config["assign_locus_numbers"] or config["assign_locus_numbers"]["completed_wizard"] != "True":

		#############################################
		# Part 4: Assign locus numbers
		#############################################
		
		print("\nThe next step is a quick one...step 4: assign_locus_numbers.\n")

		config["assign_locus_numbers"] = {}

		build = "hg38" if get_yes_no("Just one question. Are you using the hg38 build? (yes/no)\n\n", config) else "hg19"

		if build != "hg38":
			print("OK, in that case we'll assume you're using hg19. If you're using a different build or you want to use a different locus partitioning than the fault, please check out the README file for details on how to do that.\n")

		print("\nWe're already finished with the assign_locus_numbers tool.\n")

		config["assign_locus_numbers"]["genome_file"] = build
		
		config["assign_locus_numbers"]["completed_wizard"] = "True"
		
		save_results(config)	

	else:	
		print("Step 4: assign_locus_numbers has already been configured...")

	if "classify_results" not in config or "completed_wizard" not in config["classify_results"] or config["classify_results"]["completed_wizard"] != "True":

		#############################################
		# Part 5: Classify results
		#############################################

		config["classify_results"] = {}
		config["classify_results"]["rules"] = {}

		print("Now for the good stuff...part 5: classify_results.\n\nHere you decide how you want to slice and dice your results for your eventual heatmap(s).\n")

		yesno = get_yes_no("One thing we can do is sort the loci by the number of candidate QTL genes, and the number remaining after colocalization. Do you want to make such a rule? (yes/no)\n\n", config)

		if yesno:
			print("\nAwesome, good choice. I'm going to suggest using the following rules, take a look and see if this fits your needs:\n")

			pp.pprint(default_num_colocs_rule)

			keep = get_yes_no("\nDo you want to add this classification rule for your GWAS loci? (yes/no)\n\n", config)

			if not keep:
				print("\nOK, no problem. In that case, if you still want to make a custom rule of this sort, please take a look at the README section on 'num_colocs' rules, which will show you the options for defining this kind of custom rule.\n")

			else:
				rule_name = screen_input("\nWhat do you want to call your rule? (If you don't know, you could just go with 'num_colocs'). This will be the name of the new column that indicates what type of locus each test falls into.\n\n", config)
				config["classify_results"]["rules"][rule_name] = default_num_colocs_rule

		yesno = get_yes_no("\nWe can also make 'specificity' rules to classify results based on the presence or absence of colocalizations in specific tissues, traits, QTL types, or other criteria. Do you want to make such a rule? (yes/no)\n\n", config)

		if yesno:

			extended_header = header[:]
			if "mutations" in config["mutate_columns"]:
				for new_column in config["mutate_columns"]["mutations"]:
					extended_header.append(new_column["out"])

			# For each rule
			while True:

				column = screen_input(f"\nAs a reminder, your columns are {extended_header}. Which of these columns would you like to use to define your rule? (If you need information from more than one column, consider defining multiple separate rules, one at a time.)\n\n", config)

				if column not in extended_header:
					print("Sorry, that column doesn't seem to be in the header. Try again...\n\n")
					continue

				rule_name = screen_input("\nWhat do you want the name of your rule to be? (e.g. 'tissue-specificity')\n\n", config)

				config["classify_results"]["rules"][rule_name] = {"type": "specificity", "column": column, "categories": []}

				if column in header:
					with open(coloc_file) as f:
						column_vals = set([])
						f.readline()
						for line in f:
							data = line.strip().split()
							column_vals.add(data[header.index(column)])
				else:
					# This column is only being created after mutations
					mutated_column = [mc for mc in config["mutate_columns"]["mutations"] if mc["out"] == column][0]
					column_vals = set(mutated_column["map"].values())

				while True:

					cat_name = screen_input("\nWhat do you want the name of a category within this rule to be? (e.g. 'adipose-specific') Note that categories will be prioritized in order of entry, if loci meet criteria for more than one category.\n\n", config)

					category = {"class_name": cat_name}

					rule_types = ["contains_exactly", "contains_all", "contains_some", "contains_none", "contains_only"]
					
					while True:
						# Define rules one at a time, based on rules in that column...
						print(f"\nThere are a few types of criteria you can define: {rule_types}. For more information about these criteria, please view the README.\n")

						crit_type = screen_input("Which criterion type would you like to define next? (You can define more than one criterion for a single category.)\n\n", config)

						if crit_type not in rule_types:
							print("\nThat isn't one of the valid types of criteria.\n")
							continue

						category[crit_type] = []

						addition = print(f"All right, now which values do you want this {crit_type} criterion to apply to? The possible values in this column are {sorted(list(column_vals))}\n", config)
						
						while True:

							addition = screen_input("Enter one value:\n\n", config)

							if addition not in column_vals:
								print("That isn't one of the valid column values.")
								continue

							category[crit_type].append(addition)

							print("\nThe category definition right now:\n")
							pp.pprint(category)

							yesno = get_yes_no("\nDo you want to add any other values to the list for this criterion? (yes/no)\n\n", config)
							if not yesno:
								break															

						yesno = get_yes_no("\nDo you want to add any other criteria for this category? (yes/no)\n\n", config)
						if not yesno:

							config["classify_results"]["rules"][rule_name]["categories"].append(category)
							break

					print("\nGreat. Your entire rule looks like this right now:\n")
					pp.pprint(config["classify_results"]["rules"][rule_name])
					
					yesno = get_yes_no("\nDo you want to add any other categories in this rule? (yes/no)\n\n", config)
					if not yesno:

						reject = get_yes_no("\nOK, one more thing. If you've made a mistake on this rule, we can start over. Do you want to drop this rule from the config? (yes/no)\n\n", config)

						if reject:
							del config["classify_results"]["rules"][rule_name]
							print("\nRule dropped.\n")
						break

				another_rule = get_yes_no("Now, do you want to add additional rule(s)?", config)
						
				if another_rule:
					continue

				break	

		config["classify_results"]["completed_wizard"] = "True"

		print("We're now done configuring the classify_results tool.")
		
		save_results(config)	
	
	else:	
		print("Step 5: classify_results has already been configured...")

	if "make_heatmaps" not in config or "completed_wizard" not in config["make_heatmaps"] or config["make_heatmaps"]["completed_wizard"] != "True":
	
		#############################################
		# Part 6: Make heatmaps
		#############################################

		print("And now for the final tool...part 6: make_heatmaps.")

		config["make_heatmaps"] = {}

		while True:
			gwas_column = screen_input(f"One more time, the columns in your data at this point are {header}. Which column would you like to use to indicate your GWAS file? (if further customization is needed, please see the README file)", config)

			if gwas_column not in header:
				print("That column's not one of the options, try again.")
				continue

			config["make_heatmaps"]["gwas_column"] = gwas_column

			break

		while True:
			tissue_column = screen_input(f"Which column indicates your QTL tissue (or other QTL context)?", config)

			if tissue_column not in header:
				print("That column's not one of the options, try again.")
				continue

			config["make_heatmaps"]["tissue_column"] = tissue_column

			break

		yesno = get_yes_no("Do you want to specify a QTL type column? (This column should contain only the values 'eqtl' and 'sqtl'). (yes/no)", config)

		if yesno:	
			while True:
				qtl_column = screen_input(f"Which column is your QTL type column?", config)

				if gwas_column not in header:
					print("That column's not one of the options, try again.")

				config["make_heatmaps"]["qtl_column"] = qtl_column

				break

		scores_in_cells = get_yes_no("Optionally, the heatmap can display the numerical scores for each colocalization test within the cells. Do you want to do this? (yes/no)", config)

		config["make_heatmaps"]["put_scores_in_cells"] = str(scores_in_cells)

		print("Once the heatmap has been generated, you may wish to review results by collapsing across traits, tissues, genes, or some combination of all of of these. If you decide to do this, view the README file and look at the sections on the 'x_axis_collapse' and 'y_axis_collapse' parameters.")
					
		print("It's possible to make multiple sets of heatmaps using the same results, each one sorting the loci a bit differently. In any case, we'll define these sets one at a time here.")

		config["make_heatmaps"]["file_strata"] = []

		while True:

			this_stratum = {}

			out_dir = screen_input("What subdirectory would you like your heatmaps placed in, relative to the main output directory? (If you want them placed in the top level of the make_heatmaps directory, just enter nothing):", config)

			this_stratum["out_dir"] = out_dir

			concise = get_yes_no("Would you like to omit results for genes that have no colocalizations in any tested condition?", config)
			this_stratum["concise"] = str(concise)

			has_split_factor = get_yet_no("Do you want to split your heatmaps into sets based on one or more of the columns?")

			if has_split_factor:
				split_factors = []
				while True:
					factor = screen_input(f"Possible columns for splitting are {header}. Enter the name of a column you'd like to split on:", config)
					if factor not in header:
						print("Oops, I don't see that column in the header. Try again.")
						continue

					split_factors.append(factor)


					more_factors = get_yes_no(f"Currently you're splitting this set of heatmaps based on {split_factors}. Do you want to add another column to split on for this set?", config)

					if more_factors:
						continue
					break

			else:
				split_factors = ["no-split"]
			
			this_stratum["split_factors"] = split_factors
			file_strata.append(this_stratum)

			print("Great, we've added this heatmap set configuration to the config set. (If you want to modify further settings, such as excluding certain sets of results from the heatmap, see the README for info on how to make these further modifications.)")

			print("Your current heatmap configurations:")

			pp.pprint(config["make_heatmaps"])
			
			another_set = get_yes_no("Do you want to define another way of stratifying the sets of heatmaps? (yes/no)", config)

			if another_set:
				continue

			break
		
		save_results(config)	

	else:	
		print("Step 6: make_heatmaps has already been configured...")

	print("\nThat's it for the config wizard! If you want to change further details not included in this wizard, take a look at the README. Otherwise, you're all set to run the Snakemake pipeline.\n\n")

	if config_file.replace("config/", "").replace(".config", "") != coloc_file.replace("data/coloc_results/", "").replace("_colocalization_results.txt", ""):
		print("ALERT: For the pipeline to run, don't forget to rename your config or colocalization results file so that both have the same prefix (e.g. 'ir_colocalization_results.txt' and 'ir.config').")

main()

