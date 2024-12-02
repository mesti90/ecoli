#!/usr/bin/python3
'''Automatically downloading all genomes belonging to a certain species from NCBI refseq and genbank databases
'''

import csv
import re
import os
import argparse
import multiprocessing as mp
import time
from collections import defaultdict
import pandas as pd
from functools import reduce
import requests
from urllib.parse import urlparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import shutil
import logging
from itertools import chain

#Utility functions

def get_pargs():
	"""Parse command-line arguments for the script."""
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-n",help="Number of threads",default=80,type=int)
	parser.add_argument("-gd",help="Directory of genomes",default="Genomes")
	return parser.parse_args()


def init_logging():
	"""
	Initialize logging for the application.
	"""
	log_filename = "{}.{}.log".format(__file__.replace(".py3", ""), time.strftime("%Y%m%d%H%M%S"))
	logging.basicConfig(
		format='%(asctime)s. %(levelname)s: %(message)s',
		filename=log_filename,
		level=logging.INFO,
		datefmt='%Y-%m-%d %H:%M:%S'
	)
	logging.info("Run started")

def msg(text):
	"""Print timestamped messages."""
	uzenet = f"{time.strftime('%H:%M:%S')} {text}"
	print(uzenet)
	logging.info(uzenet)

def mkdir_force(path):
	"""Create a directory if it doesn't exist."""
	os.makedirs(path, exist_ok=True)


def emptyfile(out,header=None):
	if header:
		with open(out,"w") as g:
			g.write("\t".join(map(str,header))+"\n")
	else:
		if os.path.isfile(out):
			os.remove(out)
		os.system("touch " + out)



def sysexec(cmd):
	"""Execute a system command."""
	msg(f"Executing: {cmd}")
	os.system(cmd)

def read_file_lines(fname, split=False):
	"""Read lines from a file, skipping comments and duplicates."""
	msg(f"Reading {fname}")
	if not os.path.isfile(fname):
		return []
	with open(fname) as f:
		lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]
	return list(set(line.split()[0] if split else line for line in lines))

def download(link):
	"""Download a file via a given link."""
	result = link.split()[-1].rsplit(".", 1)[0]
	if os.path.isfile(result):
		return
	sysexec(link)

def wget(link):
	sysexec(f"wget {link}")


def smartint(x):
	"""	Convert a value to an integer, returning a default value for invalid inputs."""
	try:
		return int(x)
	except ValueError:
		return 0
	except:
		return -1

def load_or_create_feather(input_file, feather_file, sep="\t"):
	if os.path.exists(feather_file):
		msg(f"Loading DataFrame from Feather: {feather_file}")
		df = pd.read_feather(feather_file)
	else:
		msg(f"Reading DataFrame from CSV: {input_file}")
		df = pd.read_csv(input_file, sep=sep, low_memory=False, engine="c", dtype=str)
		msg(f"Saving DataFrame to Feather: {feather_file}")
		df.to_feather(feather_file)
	return df

#Main functions

def define_config():
	"""Define and return configuration dictionary."""
	metadata_dir = "Metadata/"
	bases = { key: f"{metadata_dir}{key}." for key in ["refseq","genbank","atb","pathogens"] }
	# Subcategories
	subcategories = {
		"ecoli": "ecoli.tsv",
		"biosample": "biosample.tsv",
		"assembly_summary": "assembly_summary.tsv"
	}
	# Configuration dictionary
	config = {
		**{f"{key}": bases[key] + subcategories["assembly_summary"] for key in bases},
		**{f"{key}_{sub}": bases[key] + subcategories[sub] for key in bases for sub in ["ecoli", "biosample"]},
		"databases": ["pathogens", "refseq", "genbank", "atb"],
		"species": "Escherichia coli",
		"atb": "Metadata/atb.E_coli_hq_GG.tsv",
		"pathogens": "Metadata/pathogens.E.coli_and_Shigella_NCBI_pathogenes.tsv",
		"biosamples": "Metadata/biosamples.all_databases.tsv",
		"biosample_metadata": "Metadata/biosamples.metadata_from_ncbi.txt",
		"biosample_metadata_tsv": "Metadata/biosamples.metadata_from_ncbi.tsv",
		"biosample_metadata_fail": "Metadata/biosamples.metadata_from_ncbi.err",
		"biosample_metadata_log": "Metadata/biosamples.metadata_from_ncbi.log",
		"pathogens_urinary": "Metadata/urinary.pathogens.tsv",
		"biosample_urinary": "Metadata/urinary.biosamples.tsv",
		"urinary": "Metadata/urinary.genomes.tsv",
		'urinary_reads': "Metadata/urinary.genomes.readlinks.tsv",
		"urinary_with_dates_and_locations": "Metadata/urinary.genomes.dates_and_locations.tsv",
		"reference": "Reference/GCF_000005845.2_ASM584v2_genomic.fna",
		"reference_link": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
		"denovo": "Denovo",
		"trimmed": "Trimmed_reads"
	}
	return config



def download_assembly_summary(source, target):
	"""Download assembly summary for the specified source."""
	if os.path.isfile(target):
		msg(f"{target} exists, skipping download.")
		return
	wget_cmd = f"wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/{source}/assembly_summary_{source}.txt -O {target}"
	sysexec(wget_cmd)

def filter_summary(source,result):
	"""Filter the assembly summary for the target species."""
	if os.path.isfile(result):
		print(f"{result} is ready, I skip filtering")
		return
	print(f"Filtering {source} -> {result}")
	with open(source) as f, open(result,"a") as g:
		rdr = csv.reader(f,delimiter='\t')
		wtr = csv.writer(g,delimiter="\t")
		for row in rdr:
			if len(row) > 7:
				block = row[7].strip("[").strip("]").split()
				if " ".join(block[:2]) == config['species']:
					wtr.writerow(row)

def get_biosamples_table(summary,output,column,header=True):
	"""Extract unique biosample IDs from the assembly summary."""
	if os.path.isfile(output):
		msg(f"Skipping biosample extraction: {output} exists.")
		return
	msg(f"Extracting biosamples from {summary} to {output}")
	with open(summary) as f:
		rdr = csv.reader(f, delimiter="\t")
		if header:
			next(rdr)
		biosamples = set(row[column] for row in rdr)
	with open(output, "w") as g:
		for bs in biosamples:
			g.write(bs + "\n")


###Currently out of use
def collect_unique_identifiers(genomes,download_file,unique_genomes_file):
	"""Collect unique genome identifiers and generate download commands."""
	if os.path.isfile(unique_genomes_file):
		msg(f"{unique_genomes_file} is already created")
		return
	msg(f"Collecting unique identifiers from {genomes}")
	with open(genomes) as f:
		rdr = csv.reader(f,delimiter="\t")
		biosamples = {}
		links = {}
		asm_name = {}
		for row in rdr:
			if row[2] not in biosamples:
				biosamples[row[2]] = []
			biosamples[row[2]].append(row[0])
			links[row[0]] = row[19]
			asm_name[row[0]] = row[15]
	unique_genomes = []
	download_cmd = []
	with open(download_file,"w") as g, open(unique_genomes_file,"w") as g2:
		wtr2 = csv.writer(g2,delimiter="\t")
		wtr2.writerow(["Genome","Biosample"])
		for sample in biosamples:
			if len(biosamples[sample]) == 1:
				best_genome = biosamples[sample][0]
			else:
				ranking = {}
				for item in biosamples[sample]:
					if item[2] == "F":
						ranking[item] = 100
					else:
						ranking[item] = 0
					ranking[item] += smartint(item[-1])
				best_genome = max(biosamples[sample],key = lambda x: ranking[x])
			link = os.path.join(f"{links[best_genome]}",f"{os.path.basename(links[best_genome])}_genomic.fna.gz")
			out_genome = os.path.join(GENOME_DIR, best_genome)
			g.write(f"wget -q -nc {link} -O {out_genome}.fna.gz\n")
			wtr2.writerow([best_genome,sample])
	msg(f"Ready with {download_file} and {unique_genomes_file}")




#Currently unused
def download_all(dl_file,downloaded):
	if os.path.isfile(downloaded):
		print("Genomes are already downloaded, skipping")
		return
	pool = mp.Pool(pargs.n)
	links = read_file_lines(dl_file)
	pool.map(download,links)
	sysexec(f"touch {downloaded}")


def merge_biosamples(databases, output):
	"""Merge biosample data from multiple databases."""
	if os.path.isfile(output):
		msg(f"{output} exists, skipping merge.")
		return
	biosample_data = {}
	all_biosamples = set()
	
	# Aggregate biosample data
	for db in databases:
		with open(config[f"{db}_biosample"]) as f:
			these_biosamples = set(line.strip() for line in f)
		biosample_data[db] = these_biosamples
		all_biosamples |= these_biosamples
	
	# Write merged biosample table
	with open(output,"w") as g:
		wtr = csv.DictWriter(g,delimiter="\t", fieldnames = ["BioSample"] + databases + ["selected_database"])
		wtr.writeheader()
		for bs in sorted(all_biosamples):
			row = {db: (1 if bs in biosample_data[db] else 0)  for db in databases}
			row["BioSample"] = bs
			wtr.writerow(row)
	print(f"Ready with {output}")




def get_urinary_pathogens(fname, output_file):
	"""
	Filter a tab-delimited file based on pathogen criteria and move the BioSample column to the first column.
	
	:param fname: Input file name (tab-delimited).
	:param output_file: Output file name for filtered results.
	"""
	if os.path.isfile(output_file):
		msg(f"I refuse to overwrite {output_file}, skipping filtering of pathogens")
		return
	msg(f"Get urinary related genomes from {fname}")
	try:
		with open(fname, "r", newline="", encoding="utf-8") as infile, open(output_file, "w", newline="", encoding="utf-8") as outfile:
			reader = csv.reader(infile, delimiter="\t")
			writer = csv.writer(outfile, delimiter="\t")
			
			# Read the header and find the BioSample column index
			header = next(reader)
			biosample_index = header.index("BioSample") if "BioSample" in header else None
			
			if biosample_index is None:
				raise ValueError("BioSample column not found in the input file header.")
			
			# Move the BioSample column to the first position in the header
			reordered_header = [header[biosample_index]] + header[:biosample_index] + header[biosample_index+1:]
			writer.writerow(reordered_header)
			
			# Filter rows based on conditions and reorder columns
			for row in reader:
				if (
					row[31] != "" and  # Column 32 is not empty
					re.search(r"uri(i?)ne|urinary|uret(h?)ra|urikult|urogenital|bladder|kidney", row[7], re.IGNORECASE) and  # Column 8 matches urinary-related terms
					not re.search(r"can(i|u|iu)s|canine|felis|feline|llama|dolphin|beef|bovine", row[7], re.IGNORECASE) and  # Column 8 does not match animal-related terms
					re.match(r"^homo|^hu|^\"homo", row[20], re.IGNORECASE) and  # Column 21 matches human-related terms
					row[48] == "562"  # Column 49 equals 562
				):
					# Move the BioSample column to the first position in the row
					reordered_row = [row[biosample_index]] + row[:biosample_index] + row[biosample_index+1:]
					writer.writerow(reordered_row)
		
		print(f"Filtered results with reordered columns saved to {output_file}")
	except Exception as e:
		print(f"Error processing file: {e}")


def fetch_biosample_data(biosample_table,output,failed,log):
	"""Fetch biosample data from NCBI and log the results."""
	if os.path.isfile(output):
		msg(f"I don't want to overwrite {output}, so I skip fetching biosample data from NCBI")
		return
	sysexec(f'''tail -n +2 {biosample_table} | cut -f 1 | xargs -P 8 -I {{}} sh -c 'efetch -db biosample -format full -id {{}} || echo "Failed ID: {{}}" >> {failed}' > {output} 2> {log}''')


def convert_raw_biosamples_to_csv(fname, output_file):
	"""Convert raw biosample data to a structured CSV file."""
	if os.path.isfile(output_file):
		msg(f"{output_file} exists, I refuse to overwrite, skipping converting raw biosample data to table")
		return
	msg(f"Convert {fname} to table")
	records = []
	my_entry = defaultdict(str)
	with open(fname) as f:
		for line in f:
			line = line.strip()
			if line.startswith("1: "):
				if my_entry:
					records.append(my_entry)
				my_entry = defaultdict(str)
				my_entry["Entry notes"] = line[3:].strip()
			elif line.startswith("Organism: "):
				my_entry['Organism'] = line[len("Organism: "):].strip()
			elif line.startswith("Identifiers: "):
				id_list = line[len("Identifiers: "):].split(";")
				for item in id_list:
					if ":" in item:
						key, value = item.strip().split(":", 1)
						my_entry[key.strip()] = value.strip()
			elif line.startswith("/") and "=" in line:
				key,value = line[1:].split("=",1)
				my_entry[key.strip("/")] = value.strip('"')
		if my_entry:
			records.append(my_entry)
			
	# Get all unique attribute keys
	all_attributes = set()
	for record in records:
		all_attributes.update(record.keys())
	
	first_cols = ['BioSample', 'Sample name', 'Organism']
	
	remaining_attributes =  all_attributes - set(first_cols)
	
	# Define column headers
	header = first_cols + sorted(remaining_attributes)
	
	# Write to CSV
	with open(output_file, 'w', newline='', encoding="utf-8") as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=header, delimiter="\t")
		writer.writeheader()
		for record in records:
			row = {attribute: record.get(attribute,'') for attribute in header}
			writer.writerow(row)
	msg(f"Ready with {output_file}")



def get_urinary_biosamples(input_file, output_file):
	"""Filter urinary-related biosamples based on specific conditions and save to a file."""
	if os.path.isfile(output_file):
		msg(f"Refusing to overwrite existing {output_file}, skipping filtering urinary genomes from biosample metadata")
		return
		
	msg(f"Get urinary related genomes from {input_file}")
	
	#Columns to prioritize in the output
	columns_to_select = ["BioSample", "collection date", "geographic location",  "isolation source", "strain", "source type", "host", "host description", "host disease", "host disease outcome", "study disease", "serotype", "pathotype", "environmental sample", "sample type", "Entry notes"]
	
	#Keywords for filtering"
	keywords = [x.lower() for x in ["CAUTI", "periurethral", "uroepithelium", "uroepithelial", "renal", "tubular epithelial", "umbrella cells", "superficial facet cells", "uroplakin", "uropathogenic", "ureter", "urination", "prostatitis", "pyelonephritis", "pyrexia", "Genitourinary", "Urinary", "Urine", "UTI", "Bladder", "Cystitis", "Catheter", "UPEC", "kidney", "urethra", "urethral", "nephron", "Nephritis"]]
	human_keywords = [x.lower() for x in ["Homo", "homo", "human"]]
	
	df = load_or_create_feather(input_file, f"{input_file}.feather")
	
	msg("Dropping rows with non-E.coli, with missing collection date or geographic location")
	conditions = [
		df['Organism'].str.startswith("Escherichia coli"),
		df['collection date'].str.strip().ne(""),
		df['collection date'].str.strip().ne("missing"),
		df['geographic location'].str.strip().ne(""),
		df['geographic location'].str.strip().ne("missing")
	]
	combined_condition = reduce(lambda x, y: x & y, conditions)
	df_filtered = df[combined_condition].copy()
	
	# Reorder columns: prioritized columns first, then the rest	
	remaining_columns = [col for col in df_filtered.columns if col not in columns_to_select]
	df_filtered = df_filtered[columns_to_select + remaining_columns].copy()
	
	
	# Search for keywords in any of the columns in `columns_to_select`
	msg("Checking urinary keywords")
	df_filtered['keyword_match'] = df_filtered[columns_to_select].apply(
	lambda row: any(keyword in " ".join(map(str, row)).lower() for keyword in keywords), axis=1
	)

	# Search for human_keywords in the entire row
	msg("Checking human keywords")
	df_filtered['human_match'] = df_filtered[columns_to_select].apply(
		lambda row: any(human_keyword in " ".join(map(str, row)).lower() for value in row for human_keyword in human_keywords),
		axis=1
	)

	# Keep rows where both keyword and human matches are true
	msg("Filtering")
	final_df = df_filtered[df_filtered['keyword_match'] & df_filtered['human_match']].drop(
		columns=['keyword_match', 'human_match']
	)
	
	final_df = final_df.dropna(axis=1, how ='all').copy()

	# Write the filtered results to a file
	final_df.to_csv(output_file, index=False, sep="\t")
	msg(f"Filtered results written to {output_file}")

def merge_pathogens_and_biosample_metadata(pathogens_file, biosamples_file, biosample_sources, output_file):
	"""Merge tables by the 'BioSample' column and add prefixes to distinguish column origins."""
	if os.path.isfile(output_file):
		msg(f"{output_file} is ready, skip merging")
		return
	msg(f"Merging {pathogens_file} and {biosamples_file}")
	# Read the tables
	pathogens_df = pd.read_csv(pathogens_file, sep="\t", dtype=str)
	biosamples_df = pd.read_csv(biosamples_file, sep="\t", dtype=str)
	sources_df = pd.read_csv(biosample_sources, sep="\t", dtype=str)
	

	# Add prefixes to column names (except 'BioSample')
	pathogens_df = pathogens_df.add_prefix("Pathogens|")
	biosamples_df = biosamples_df.add_prefix("Biosamples|")

	# Restore 'BioSample' column (it gets prefixed by add_prefix)
	pathogens_df = pathogens_df.rename(columns={"Pathogens|BioSample": "BioSample"})
	biosamples_df = biosamples_df.rename(columns={"Biosamples|BioSample": "BioSample"})

	# Merge the tables
	merged_df = pd.merge(
		biosamples_df,
		pathogens_df,
		on="BioSample",
		how="outer"  # Use 'outer' to keep all rows
	)
	merged_df = pd.merge(sources_df, merged_df, on = "BioSample", how= "right")

	# Save the merged DataFrame
	merged_df.to_csv(output_file, sep="\t", index=False)
	msg(f"Merged table saved to {output_file}")


def fetch_read_identifiers(infile, batch_size=20):
	"""Download genome information in batches of BioSamples."""
	# Read the BioSample IDs from the input file
	mkdir_force("Metadata/Reads")
	downloaded_flag = "Metadata/Reads/batches_downloaded.flag"
	if os.path.isfile(downloaded_flag):
		msg("Read identifiers are already retrieved. Skipping this part.")
		return
	biosample_ids = get_biosamples(infile)
	# Create batches of BioSamples
	for i in range(0, len(biosample_ids), batch_size):
		# Create a query with up to `batch_size` BioSamples
		batch = biosample_ids[i:i + batch_size]
		output_file = f"Metadata/Reads/genomes_batch_{i//batch_size + 1}.tsv"
		if os.path.isfile(output_file):
			continue
		query = " OR ".join([f"sample_accession={sample}" for sample in batch])

		# Construct the API URL
		api_url = f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=({query})&fields=sample_accession,instrument_platform,fastq_ftp"

		# Perform the request
		msg(f"Fetching data? API url: {api_url}")
		response = requests.get(api_url)

		# Check for successful response
		if response.status_code == 200:
			with open(output_file, "w") as f:
				f.write(response.text)
			print(f"Saved batch {i//batch_size + 1} to {output_file}")
		else:
			print(f"Failed to fetch batch {i//batch_size + 1}: {response.status_code}")
	emptyfile(downloaded_flag)


def get_biosamples(infile):
	return pd.read_csv(infile, sep="\t", usecols=["BioSample"])["BioSample"].tolist()

def concat_batches(infile, output_file, batch_size=20):
	"""Concatenate all batch files into a single file, keeping only one header."""
	if os.path.isfile(output_file):
		msg(f"{output_file} exists. I'll skip concatenating batches.")
		return
	biosample_ids = get_biosamples(infile)
	# Create batches of BioSamples
	temp_files = [f"Metadata/Reads/genomes_batch_{i//batch_size + 1}.tsv" for i in range(0, len(biosample_ids), batch_size)]

	dfs = []
	header = []
	for i, temp_file in enumerate(temp_files):
		if not os.path.isfile(temp_file):
			continue
		df = pd.read_csv(temp_file, sep="\t")
		if not header:
			df = df.iloc[1:]  # Skip the header for subsequent files
		dfs.append(df)

	# Combine all DataFrames
	merged_df = pd.concat(dfs, ignore_index=True)

	# Save the final concatenated DataFrame
	merged_df.to_csv(output_file, sep="\t", index=False)
	print(f"All batch files concatenated into {output_file}")



def download_reads(input_file, genomes, output_dir="Raw_reads", max_threads=20, errors="Metadata/biosamples.no_reads_from_ena.list"):
	"""
	Download reads from the provided file and rename them based on BioSample IDs.
	
	:param input_file: Path to the file containing run_accession, sample_accession, and fastq_ftp.
	:param output_dir: Directory where the downloaded reads will be stored.
	"""
	# Ensure the output directory exists
	flag = "Metadata/raw_reads_downloaded.flag"
	#if os.path.isfile(flag):
	#	msg("All reads are downloaded")
	#	return
	msg("Downloading raw reads - collecting download links")
	os.makedirs(output_dir, exist_ok=True)
	biosample_list = get_biosamples(genomes)
	download_tasks = []
	df = pd.read_csv(input_file, sep="\t").fillna("")
	not_illumina = 0
	runs = 0
	for _,row in df.iterrows():
		if row['sample_accession'] not in biosample_list:
			logging.info(f"Skipping {row['sample_accession']} as it's not among the selected genomes")
			continue
		if row["instrument_platform"] != "ILLUMINA":
			logging.info(f"Skipping {row['sample_accession']} as it's not sequenced on Illumina platform")
			not_illumina += 1
			continue
		if row['fastq_ftp']:
			runs += 1
			for item in row['fastq_ftp'].split(";"):
				download_tasks.append([item, os.path.join(output_dir, os.path.basename(item))])
	msg(f"Now starting download for {len(download_tasks)} links for {runs} runs. Skipping {not_illumina} runs which are not sequenced on Illumina")
	with ThreadPoolExecutor(max_threads) as executor:
		executor.map(lambda task: download_file(*task), download_tasks)
	emptyfile(flag)



def download_file(link, output_path):
	"""
	Download a file from the given link and save it to the specified path.
	
	:param link: URL of the file to download.
	:param output_path: Path to save the downloaded file.
	"""
	temp_path = f"{output_path}.part"
	if os.path.isfile(output_path) or os.path.isfile(temp_path):
		msg(f"{output_path} is already downloaded or being downloaded")
		return
	if not link.startswith("http"):
		link = f"https://{link}"
	msg(f"Downloading {link} to {output_path}...")
	try:
		response = requests.get(link, stream=True)
		response.raise_for_status()
		with open(temp_path, "wb") as temp_file:
			for chunk in response.iter_content(chunk_size=8192):
				temp_file.write(chunk)
		os.rename(temp_path, output_path)
		print(f"Downloaded and saved: {output_path}")
	except requests.exceptions.RequestException as e:
		print(f"HTTP request failed for {link}: {e}")
	except Exception as e:
		print(f"An error occurred while downloading {link}: {e}")
	finally:
		try:
			if os.path.isfile(temp_path):
				os.remove(temp_path)
		except Exception as cleanup_error:
			print(f"Failed to clean up temp file {temp_path}: {cleanup_error}")



def filter_and_transform_genomes(input_file, output_file):
	"""Filter and transform the urinary_genomes.tsv file."""
	# Read the input file
	df = pd.read_csv(input_file, sep="\t", dtype=str)

	# Define invalid values for dates and locations
	invalid_dates_biosample = [
		"0000", "Missing", "not applicable", "Not applicable",
		"not collected", "not provided", "unknown", "Unknown"
	]
	invalid_dates_pathogens = ["unknown", "Unknown", ""]

	invalid_locations = ["", "unknown", "Unknown"]

	# Process Biosamples|collection date
	df["Biosamples|collection date"] = df["Biosamples|collection date"].apply(
		lambda x: x if (
			pd.notna(x)
			and x not in invalid_dates_biosample
			and not "/" in x
			and x.isdigit()
			and int(x) >= 1950
		) else None
	)

	# Process Pathogens|Collection date
	df["Pathogens|Collection date"] = df["Pathogens|Collection date"].apply(
		lambda x: x if (
			pd.notna(x)
			and x not in invalid_dates_pathogens
			and x.isdigit()
			and int(x) >= 1950
		) else None
	)

	# Create collection_date column
	df["collection_date"] = df.apply(
		lambda row: row["Biosamples|collection date"]
		if pd.notna(row["Biosamples|collection date"])
		else row["Pathogens|Collection date"]
		if pd.notna(row["Pathogens|Collection date"])
		else None,
		axis=1
	)

	# Drop rows where collection_date is None
	df = df[df["collection_date"].notna()]

	# Process Biosamples|geographic location
	df["Biosamples|geographic location"] = df["Biosamples|geographic location"].apply(
		lambda x: x if (
			pd.notna(x)
			and x not in invalid_locations
		) else None
	)

	# Process Pathogens|Location
	df["Pathogens|Location"] = df["Pathogens|Location"].apply(
		lambda x: x if (
			pd.notna(x)
			and x not in invalid_locations
		) else None
	)

	# Create geographic_location column
	df["geographic_location"] = df.apply(
		lambda row: (
			f"{row['Biosamples|geographic location']}|{row['Biosamples|geographic location (region and locality)']}"
			if pd.notna(row["Biosamples|geographic location (region and locality)"])
			else row["Biosamples|geographic location"]
			if pd.notna(row["Biosamples|geographic location"])
			else row["Pathogens|Location"]
			if pd.notna(row["Pathogens|Location"])
			else None
		),
		axis=1
	)

	# Drop rows where geographic_location is None
	df = df[df["geographic_location"].notna()]

	# Rearrange columns: Put the new columns first
	first_cols = ["BioSample","pathogens","refseq","genbank","atb","collection_date", "geographic_location"]
	column_order = first_cols + [
		col for col in df.columns if col not in first_cols
	]
	df = df[column_order]

	# Save the filtered DataFrame to the output file
	df.to_csv(output_file, sep="\t", index=False)
	print(f"Filtered and transformed data saved to {output_file}")


def download_and_extract(url, output_file):
	"""Download a gzip file and extract it on the fly."""
	try:
		response = requests.get(url, stream=True)
		response.raise_for_status()  # Raise an error for HTTP issues
		with gzip.open(response.raw, "rb") as gz_file:
			with open(output_file, "wb") as f_out:
				shutil.copyfileobj(gz_file, f_out)
		msg(f"Downloaded and extracted {url} to {output_file}")
	except requests.exceptions.RequestException as e:
		msg(f"Failed to download {url}: {e}")
	except OSError as e:
		msg(f"Failed to extract {url}: {e}")

def get_reference_genome():
	if os.path.isfile(config['reference']):
		msg(f"{config['reference']} is already downloaded")
		return
	if os.path.dirname(config['reference']):
		mkdir_force(os.path.dirname(config['reference']))
	download_and_extract(config['reference_link'], config['reference'])

def prepare_assembly_commands(genomes, links):
	"""
	Generate assembly commands using cutadapt for trimming and SPAdes for assembly.

	:param genomes: Path to the genomes file containing biosamples.
	:param links: Path to the links file containing sample_accession and fastq_ftp columns.
	:param config: Dictionary containing output paths for trimmed reads and assemblies.
	"""
	msg("Generating assembly commands")
	# Get biosamples and their corresponding fastq files
	biosamples = get_biosamples(genomes)
	df = pd.read_csv(links, sep="\t", usecols=["sample_accession", "fastq_ftp"]).fillna("")
	biosample_fastq_dict = {row["sample_accession"]: get_downloaded_names(row["sample_accession"], row["fastq_ftp"], savelinks=False) for _, row in df.iterrows()}
	
	commands = []

	for bs in biosamples:
		fq_files = biosample_fastq_dict.get(bs, [])
		if not fq_files:
			logging.info(f"{bs}: missing reads, assembly not possible")
			continue

		# Classify the fastq files
		se_reads = [fq for fq in fq_files if ".SE" in fq]
		r1_reads = [fq for fq in fq_files if ".R1." in fq]
		r2_reads = [fq for fq in fq_files if ".R2." in fq]

		# Trimming commands
		trimming_commands = []
		for idx, se_read in enumerate(se_reads, start=1):
			se_trimmed = f"{config['trimmed']}/{bs}.SE{idx}.trimmed.fastq.gz"
			cmd = f"cutadapt --trim-n --quality-base 33 --max-n 0.5 -m 30 --quality-cutoff 20 -o {se_trimmed} {se_read} --quiet -Z"
			trimming_commands.append(cmd)

		for idx, (r1_read, r2_read) in enumerate(zip(r1_reads, r2_reads), start=1):
			r1_trimmed = f"{config['trimmed']}/{bs}.PE{idx}.R1.trimmed.fastq.gz"
			r2_trimmed = f"{config['trimmed']}/{bs}.PE{idx}.R2.trimmed.fastq.gz"
			cmd = f"cutadapt --trim-n --quality-base 33 --max-n 0.5 -m 30 --quality-cutoff 20 -o {r1_trimmed} -p {r2_trimmed} {r1_read} {r2_read} --quiet -Z"
			trimming_commands.append(cmd)
		
		
		
		# Assembly command
		spades_output = f"{config['denovo']}/{bs}_refguided_denovo"
		spades_log = f"{config['denovo']}/{bs}.spades.log"
		spades_err = f"{config['denovo']}/{bs}.spades.err"
		se_trimmed_reads = " ".join([f"{config['trimmed']}/{bs}.SE{idx}.trimmed.fastq.gz" for idx in range(1, len(se_reads) + 1)])
		pe_r1_trimmed_reads = " ".join([f"{config['trimmed']}/{bs}.PE{idx}.R1.trimmed.fastq.gz" for idx in range(1, len(r1_reads) + 1)])
		pe_r2_trimmed_reads = " ".join([f"{config['trimmed']}/{bs}.PE{idx}.R2.trimmed.fastq.gz" for idx in range(1, len(r2_reads) + 1)])
		assembly_cmd = (
			f"spades.py -o {spades_output} "
			+ (f"-1 {pe_r1_trimmed_reads} " if pe_r1_trimmed_reads else "")
			+ (f"-2 {pe_r2_trimmed_reads} " if pe_r2_trimmed_reads else "")
			+ (f"-s {se_trimmed_reads} " if se_trimmed_reads else "")
			+ f"--untrusted-contigs {config['reference']} "
			+ f"> {spades_log} 2> {spades_err}"
		)
		combined_cmd = f"{' && '.join(trimming_commands)} && {assembly_cmd}"
		commands.append(combined_cmd)

	return commands

def write_commands(commands,commands_file):
	with open(commands_file, "w") as cmd_file:
		for cmd in commands:
			cmd_file.write(f"{cmd}\n")
	print(f"Assembly commands written to {commands_file}")



def main():
	"""Main execution flow"""
	global config
	pargs = get_pargs()
	config = define_config()
	
	# Step 1: Download and filter assembly summaries
	for db in ['refseq','genbank']:
		download_assembly_summary(db,config[db])
		filter_summary(config[db],config[f"{db}_ecoli"])
		get_biosamples_table(config[f"{db}_ecoli"], config[f"{db}_biosample"], 2)
		
	
	# Step 2: Process manually downloaded pathogen data
	get_biosamples_table(config["pathogens"],config["pathogens_biosample"], 13)
	get_biosamples_table(config["atb"],config["atb_biosample"], 0)
	
	# Step 3: Merge biosample data
	merge_biosamples(config["databases"], config["biosamples"])
	
	
	# Step 4: Fetch, convert, and filter biosample metadata

	fetch_biosample_data(config["biosamples"], config["biosample_metadata"], config["biosample_metadata_fail"], config["biosample_metadata_log"])
	get_urinary_pathogens(config["pathogens"], config["pathogens_urinary"])
	convert_raw_biosamples_to_csv(config["biosample_metadata"], config["biosample_metadata_tsv"])
	get_urinary_biosamples(config["biosample_metadata_tsv"], config["biosample_urinary"])
	
	# Step 5: Merge tables and collect all metadata
	merge_pathogens_and_biosample_metadata(config['pathogens_urinary'], config['biosample_urinary'], config['biosamples'], config['urinary'])
	
	filter_and_transform_genomes(config['urinary'], config["urinary_with_dates_and_locations"])
	
	
	
	fetch_read_identifiers(config["urinary_with_dates_and_locations"])
	concat_batches(config["urinary_with_dates_and_locations"], config['urinary_reads'])
	#
	download_reads(config['urinary_reads'], config["urinary_with_dates_and_locations"])
	#
	#!get_reference_genome()
	#
	#commands = prepare_assembly_commands(config["urinary_with_dates_and_locations"], config['urinary_reads'])
	#write_commands(commands,"assembly_commands.sh")




	
if __name__ == "__main__":
	init_logging()
	main()
