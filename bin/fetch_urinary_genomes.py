#!/usr/bin/python3
'''Automatically downloading all genomes belonging to a certain species from NCBI refseq and genbank databases
'''
#docker: klebsiella

#TODO: download should be moved to the process_sample part

import csv
import glob
import re
import os
import argparse
import multiprocessing as mp
import time
from collections import defaultdict, OrderedDict
import pandas as pd
from functools import reduce
import requests
from urllib.parse import urlparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import shutil
import logging
from itertools import chain
import subprocess
from pathlib import Path

#Utility functions

def get_pargs():
	"""Parse command-line arguments for the script."""
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-n",help="Number of threads",default=20,type=int)
	parser.add_argument("-gd",help="Directory of genomes",default="Genomes")
	return parser.parse_args()


def init_logging():
	"""
	Initialize logging for the application.
	"""
	log_dir = "Logging"
	os.makedirs(log_dir, exist_ok=True)
	log_filename = os.path.join(log_dir, os.path.basename(f"{__file__.replace('.py3', '')}.{time.strftime('%Y%m%d%H%M%S')}.log"))
	logging.basicConfig(
		format='%(asctime)s. %(levelname)s: %(message)s',
		filename=log_filename,
		level=logging.INFO,
		datefmt='%Y-%m-%d %H:%M:%S'
	)
	print(f"Run started, logging to {log_filename}")
	logging.info(f"Run started")

def msg(text,onscreen=False):
	"""Print timestamped messages."""
	uzenet = f"{time.strftime('%H:%M:%S')} {text}"
	if onscreen:
		print(uzenet)
	logging.info(uzenet)

def info(text):
	msg(text, True)

def error(text,onscreen=True):
	"""Print timestamped messages."""
	uzenet = f"{time.strftime('%H:%M:%S')} {text}"
	if onscreen:
		print(uzenet)
	logging.error(uzenet)


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
		"reference": "Reference/GCF_000005845.2.fa",
		"reference_link": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
		"denovo": "Intermediate_files/Denovo",
		"fasterq": "Intermediate_files/Fasterq_dump",
		"assemblies": "Assemblies",
		"stats": "Statistics",
		"raw_reads": "Raw_reads",
		"trimmed": "Intermediate_files/Trimmed_reads",
		"commands": "Metadata/assembly_commands.tsv",
		'minlen': 1000,
		'mincov': 10,
	}
	return config



def download_assembly_summary(source, target):
	"""Download assembly summary for the specified source."""
	if os.path.isfile(target):
		info(f"{target} exists, skipping download.")
		return
	wget_cmd = f"wget -nc https://ftp.ncbi.nlm.nih.gov/genomes/{source}/assembly_summary_{source}.txt -O {target}"
	sysexec(wget_cmd)

def filter_summary(source,result):
	"""Filter the assembly summary for the target species."""
	if os.path.isfile(result):
		info(f"{result} is ready, I skip filtering")
		return
	info(f"Filtering {source} -> {result}")
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
		info(f"Skipping biosample extraction: {output} exists.")
		return
	info(f"Extracting biosamples from {summary} to {output}")
	with open(summary) as f:
		rdr = csv.reader(f, delimiter="\t")
		if header:
			next(rdr)
		biosamples = set(row[column] for row in rdr)
	with open(output, "w") as g:
		for bs in biosamples:
			g.write(bs + "\n")


def merge_biosamples(databases, output):
	"""Merge biosample data from multiple databases."""
	if os.path.isfile(output):
		info(f"{output} exists, skipping merge.")
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
	info(f"Ready with {output}")




def get_urinary_pathogens(fname, output_file):
	"""
	Filter a tab-delimited file based on pathogen criteria and move the BioSample column to the first column.
	
	:param fname: Input file name (tab-delimited).
	:param output_file: Output file name for filtered results.
	"""
	if os.path.isfile(output_file):
		info(f"I refuse to overwrite {output_file}, skipping filtering of pathogens")
		return
	info(f"Get urinary related genomes from {fname}")
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
		info(f"I don't want to overwrite {output}, so I skip fetching biosample data from NCBI")
		return
	sysexec(f'''tail -n +2 {biosample_table} | cut -f 1 | xargs -P 8 -I {{}} sh -c 'efetch -db biosample -format full -id {{}} || echo "Failed ID: {{}}" >> {failed}' > {output} 2> {log}''')


def convert_raw_biosamples_to_csv(fname, output_file):
	"""Convert raw biosample data to a structured CSV file."""
	if os.path.isfile(output_file):
		info(f"{output_file} exists, skipping converting raw biosample data to table")
		return
	info(f"Convert {fname} to table")
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
	info(f"Ready with {output_file}")



def get_urinary_biosamples(input_file, output_file):
	"""Filter urinary-related biosamples based on specific conditions and save to a file."""
	if os.path.isfile(output_file):
		info(f"Refusing to overwrite existing {output_file}, skipping filtering urinary genomes")
		return
		
	info(f"Get urinary related genomes from {input_file}")
	
	#Columns to prioritize in the output
	columns_to_select = ["BioSample", "collection date", "geographic location",  "isolation source", "strain", "source type", "host", "host description", "host disease", "host disease outcome", "study disease", "serotype", "pathotype", "environmental sample", "sample type", "Entry notes"]
	
	#Keywords for filtering"
	keywords = [x.lower() for x in ["CAUTI", "periurethral", "uroepithelium", "uroepithelial", "renal", "tubular epithelial", "umbrella cells", "superficial facet cells", "uroplakin", "uropathogenic", "ureter", "urination", "prostatitis", "pyelonephritis", "pyrexia", "Genitourinary", "Urinary", "Urine", "UTI", "Bladder", "Cystitis", "Catheter", "UPEC", "kidney", "urethra", "urethral", "nephron", "Nephritis"]]
	human_keywords = [x.lower() for x in ["Homo", "homo", "human"]]
	
	df = load_or_create_feather(input_file, f"{input_file}.feather")
	
	info("Dropping rows with non-E.coli, with missing collection date or geographic location")
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
	info("Checking urinary keywords")
	df_filtered['keyword_match'] = df_filtered[columns_to_select].apply(
	lambda row: any(keyword in " ".join(map(str, row)).lower() for keyword in keywords), axis=1
	)

	# Search for human_keywords in the entire row
	info("Checking human keywords")
	df_filtered['human_match'] = df_filtered[columns_to_select].apply(
		lambda row: any(human_keyword in " ".join(map(str, row)).lower() for value in row for human_keyword in human_keywords),
		axis=1
	)

	# Keep rows where both keyword and human matches are true
	info("Filtering")
	final_df = df_filtered[df_filtered['keyword_match'] & df_filtered['human_match']].drop(
		columns=['keyword_match', 'human_match']
	)
	
	final_df = final_df.dropna(axis=1, how ='all').copy()
	
	# Write the filtered results to a file
	final_df.to_csv(output_file, index=False, sep="\t")
	info(f"Filtered results written to {output_file}")

def merge_pathogens_and_biosample_metadata(pathogens_file, biosamples_file, biosample_sources, output_file):
	"""Merge tables by the 'BioSample' column and add prefixes to distinguish column origins."""
	if os.path.isfile(output_file):
		info(f"{output_file} is ready, skip merging")
		return
	info(f"Merging {pathogens_file} and {biosamples_file}")
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
	info(f"Merged table saved to {output_file}")


def fetch_read_identifiers(infile, batch_size=20):
	"""Download genome information in batches of BioSamples."""
	# Read the BioSample IDs from the input file
	mkdir_force("Metadata/Reads")
	downloaded_flag = "Metadata/Reads/batches_downloaded.flag"
	if os.path.isfile(downloaded_flag):
		info("Read identifiers are already retrieved. Skipping this part.")
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
		msg(f"Fetching data: API url: {api_url}")
		response = requests.get(api_url)

		# Check for successful response
		if response.status_code == 200:
			with open(output_file, "w") as f:
				f.write(response.text)
			msg(f"Saved batch {i//batch_size + 1} to {output_file}")
		else:
			error(f"Failed to fetch batch {i//batch_size + 1}: {response.status_code}")
	emptyfile(downloaded_flag)


def get_biosamples(infile):
	return pd.read_csv(infile, sep="\t", usecols=["BioSample"])["BioSample"].tolist()

def concat_batches(infile, output_file, batch_size=20):
	"""Concatenate all batch files into a single file, keeping only one header."""
	if os.path.isfile(output_file):
		info(f"{output_file} exists. I'll skip concatenating batches.")
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
	info(f"All batch files concatenated into {output_file}")



def download_reads(input_file, genomes, output_dir, max_threads=20, errors="Metadata/biosamples.no_reads_from_ena.list"):
	"""
	Download reads from the provided file and rename them based on BioSample IDs.
	
	:param input_file: Path to the file containing run_accession, sample_accession, and fastq_ftp.
	:param output_dir: Directory where the downloaded reads will be stored.
	"""
	# Ensure the output directory exists
	flag = "Metadata/raw_reads_downloaded.flag"
	if os.path.isfile(flag):
		info("All reads are downloaded")
		return
	info("Downloading raw reads - collecting download links")
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
	info(f"Now starting download for {len(download_tasks)} links for {runs} runs. Skipping {not_illumina} runs which are not sequenced on Illumina")
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
		time.sleep(1)
		response = requests.get(link, stream=True)
		response.raise_for_status()
		with open(temp_path, "wb") as temp_file:
			for chunk in response.iter_content(chunk_size=8192):
				temp_file.write(chunk)
		os.rename(temp_path, output_path)
		msg(f"Downloaded and saved: {output_path}")
	except requests.exceptions.RequestException as e:
		error(f"HTTP request failed for {link}: {e}")
	except Exception as e:
		error(f"An error occurred while downloading {link}: {e}")
	finally:
		try:
			if os.path.isfile(temp_path):
				os.remove(temp_path)
		except Exception as cleanup_error:
			error(f"Failed to clean up temp file {temp_path}: {cleanup_error}")



def filter_and_transform_genomes(input_file, output_file):
	"""Filter and transform the urinary_genomes.tsv file."""
	if os.path.isfile(output_file):
		info(f"Skip filtering genomes, as {output_file} exists")
		return
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
	info(f"Filtered and transformed data saved to {output_file}")


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


def check_gz_integrity(gz):
	"""Check the integrity of a gzip file."""
	if not os.path.isfile(gz):
		logging.error(f"{gz} does not exist.")
		return False
	res = os.system(f"gzip -t {gz} > /dev/null 2>&1") == 0
	msg(f"{gz} is ok" if res else f"ERROR: {gz} is corrupt")
	return res

def read_cmd(cmd):
	"""
	Execute a shell command and return its output.

	:param cmd: The shell command to execute.
	:return: The output of the command as a string, with leading and trailing whitespace removed.
	:raises RuntimeError: If the command fails to execute.
	"""
	try:
		return subprocess.check_output(cmd, shell=True, text=True).strip()
	except subprocess.CalledProcessError as e:
		raise RuntimeError(f"Command failed with error code {e.returncode}: {cmd}\n{e.output}")
	except Exception as e:
		raise RuntimeError(f"An unexpected error occurred while running the command: {cmd}\n{e}")

def get_fa_names(fa):
	"""
	Extract sequence names from a FASTA/FASTQ file using seqtk.

	:param fa: Path to the FASTA/FASTQ file.
	:return: A list of sequence names.
	:raises RuntimeError: If the command fails or the input file is invalid.
	"""
	if not os.path.isfile(fa):
		raise RuntimeError(f"File not found: {fa}")

	try:
		# Use seqtk to get sequence names
		output = read_cmd(f'seqtk comp {fa} | cut -f 1')
		return output.split("\n") if output else []
	except RuntimeError as e:
		raise RuntimeError(f"Failed to extract sequence names from {fa}: {e}")


def check_paired(fa1, fa2):
	names1 = get_fa_names(fa1)
	names2 = get_fa_names(fa2)
	return names1 == names2


def prepare_assembly_commands(genomes, links, outfile):
	"""
	Generate assembly commands using cutadapt for trimming and SPAdes for assembly.

	:param genomes: Path to the genomes file containing biosamples.
	:param links: Path to the links file containing sample_accession and fastq_ftp columns.
	:param config: Dictionary containing output paths for trimmed reads and assemblies.
	"""
	if os.path.isfile(outfile):
		info(f"{outfile} exists, I skip preparing assembly commands")
		return
	info("Generating assembly commands")
	# Get biosamples and their corresponding fastq files
	biosamples = get_biosamples(genomes)
	df = pd.read_csv(links, sep="\t").fillna("")
	
	header = ["BioSample","Run","fastq_ftp","SE", "R1", "R2","SE trim","R1 trim", "R2 trim","Same name in paired-end","trim","spades"]
	
	datarows = {}
	
	for _,row in df.iterrows():
		run_acc = row.get("run_accession")
		bs = row.get("sample_accession")
		fq_files = row.get("fastq_ftp","").split(";")
		if row.get("instrument_platform") != "ILLUMINA" or not (run_acc and bs and fq_files and bs in biosamples):
			continue
		
		# Classify the fastq files
		r1_reads = [os.path.basename(fq) for fq in fq_files if fq.endswith("_1.fastq.gz")]
		r2_reads = [os.path.basename(fq) for fq in fq_files if fq.endswith("_2.fastq.gz")]
		se_reads = [os.path.basename(fq) for fq in fq_files if fq and os.path.basename(fq) not in (r1_reads + r2_reads)]
		r1_trimmed_reads = []
		r2_trimmed_reads = []
		se_trimmed_reads = []
		
		#Map raw and trimmed file paths
		raw = {read: os.path.join(config['raw_reads'], read) for read in r1_reads + r2_reads + se_reads}
		trimmed = {read: os.path.join(config['trimmed'], f"{read.replace('.fastq.gz','')}.trimmed.fastq.gz") for read in r1_reads + r2_reads + se_reads}
		
		data = {
			"BioSample": bs,
			"Run": run_acc,
			"fastq_ftp": row.get("fastq_ftp",""),
			"SE": ", ".join(raw[read] for read in se_reads),
			"R1": ", ".join(raw[read] for read in r1_reads),
			"R2": ", ".join(raw[read] for read in r2_reads),
			"SE trim": ", ".join(trimmed[read] for read in se_reads),
			"R1 trim": ", ".join(trimmed[read] for read in r1_reads),
			"R2 trim": ", ".join(trimmed[read] for read in r2_reads),
		}
		
		# Trimming commands
		trimming_commands = [
			f"cutadapt --trim-n --quality-base 33 --max-n 0.5 -m 30 --quality-cutoff 20 "
			f"-o {trimmed[se_read]} {raw[se_read]} --quiet -Z"
			for se_read in se_reads
		] + [
			f"cutadapt --trim-n --quality-base 33 --max-n 0.5 -m 30 --quality-cutoff 20 "
			f"-o {trimmed[r1_read]} -p {trimmed[r2_read]} {raw[r1_read]} {raw[r2_read]} --quiet -Z"
			for r1_read, r2_read in zip(r1_reads, r2_reads)
		]
		
		data["trim"] = ' && '.join(trimming_commands)
		
		# Assembly command
		spades_output = os.path.join(config['denovo'], f"{bs}.{run_acc}_denovo")
		spades_log = os.path.join(config['denovo'], f"{bs}.{run_acc}.spades.log")
		spades_err = os.path.join(config['denovo'], f"{bs}.{run_acc}.spades.err")
		se_trimmed_reads = " ".join(trimmed[read] for read in se_reads)
		r1_trimmed_reads = " ".join(trimmed[read] for read in r1_reads)
		r2_trimmed_reads = " ".join(trimmed[read] for read in r2_reads)
		data["spades"] = (
			f"spades.py -t 1 -o {spades_output} "
			+ (f"-1 {r1_trimmed_reads} " if r1_reads else "")
			+ (f"-2 {r2_trimmed_reads} " if r2_reads else "")
			+ (f"-s {se_trimmed_reads} " if se_reads else "")
			+ f"--untrusted-contigs {config['reference']} > {spades_log} 2> {spades_err}"
		)
		datarows[run_acc] = data
		
	with open(outfile, "w") as g:
		wtr = csv.DictWriter(g, fieldnames=header, delimiter="\t")
		wtr.writeheader()
		
		for row in datarows.values():
			# Fill missing fields with empty strings
			complete_row = {field: row.get(field, "") for field in header}
			wtr.writerow(complete_row)
	
	msg(f"Ready with {outfile}")


def check_set_of_fq(sample, se, r1, r2):
	"""
	Check the integrity of a set of FASTQ files and ensure all required files are present.

	:param sample: Sample name for logging.
	:param se: Path to single-end FASTQ file.
	:param r1: Path to paired-end R1 FASTQ file.
	:param r2: Path to paired-end R2 FASTQ file.
	:return: True if all files are valid; False otherwise.
	"""
	files = {"SE": se, "R1": r1, "R2": r2}
	missing_files = [name for name, fq in files.items() if fq and not os.path.isfile(fq)]
	corrupt_files = [name for name, fq in files.items() if fq and os.path.isfile(fq) and not check_gz_integrity(fq)]

	# Log missing files
	if missing_files:
		error(f"{sample}: Missing files: {', '.join(missing_files)}\n")

	# Log corrupt files
	if corrupt_files:
		error(f"{sample}: Corrupt files: {', '.join(corrupt_files)}\n")

	# Return False if any issues were found
	if missing_files or corrupt_files:
		return False

	# Ensure paired-end files are correctly paired
	if (r1 and not r2) or (r2 and not r1):
		error(f"{sample}: Missing one of the paired-end files (R1 or R2)\n")
		return False

	# Ensure at least one valid file is present
	if not (se or (r1 and r2)):
		error(f"{sample}: No valid FASTQ files provided")
		return False

	return True

def collect_scaffolds(sample, assembly, denovo_path):
	"""
	Collect and filter scaffolds for a given sample based on length and coverage criteria.

	:param sample: Name of the sample.
	:param assembly: Path to the final assembly file.
	:param denovo_path: Path to the SPAdes denovo directory.
	:return: True if the assembly file is successfully created or already exists.
	"""
	msg(f"Collecting scaffolds for {sample}")
	
	# Paths for scaffolds and filtered list
	scaffolds_source = os.path.join(denovo_path,"scaffolds.fasta")
	scaffolds = f"{os.path.join(config['denovo'],sample)}_scaffolds.fasta"
	filterlist = f"{os.path.join(config['denovo'],sample)}_scaffolds.filtered.list"
	
	# Handle gzipped assembly
	gz = False
	if assembly.endswith(".gz"):
		gz = True
		assembly = assembly.rsplit(".", 1)[0]  # Remove ".gz" extension for processing
	
	# Ensure scaffolds file exists
	if not os.path.isfile(scaffolds):
		if os.path.isfile(scaffolds_source):
			shutil.copy(scaffolds_source,scaffolds)
			msg(f"{sample}: Scaffolds copied from {scaffolds_source}")
		else:
			error(f"{sample}: Scaffolds source missing at {scaffolds_source}")
			return False
	
	# Extract scaffold sequence names
	sequences = get_fa_names(scaffolds)
	if not sequences:
		error(f"{sample}: No sequences found in {scaffolds}")
		return False
	
	# Filter scaffolds based on length and coverage
	with open(filterlist,"w") as g:
		filtered_count = 0
		for seqname in sequences:
			block = seqname.split("_")
			if len(block) < 6:
				msg(f"{sample}: Invalid sequence name format: {seqname}")
				continue
			
			length = int(block[3])
			cov = float(block[5])
			if length >= config['minlen'] and cov >= config['mincov']:
				g.write(f"{seqname}\n")
				filtered_count += 1
	
	msg(f"{sample}: {filtered_count} contigs passed the filter.")
	
	# Create filtered assembly
	sysexec(f"seqtk subseq {scaffolds} {filterlist} > {assembly}")
	if gz:
		sysexec(f"gzip -f {assembly}")
		msg(f"{sample}: Gzipped assembly created at {assembly}.gz")
		
	else:
		msg(f"{sample}: assembly created at {assembly}")
	
	return True

def count_reads_and_bases(fastq):
	"""Count the number of reads and total nucleotides in a FASTQ file."""
	if not fastq or not os.path.isfile(fastq):
		return 0, 0
	cmd = f"zcat {fastq} | awk 'NR % 4 == 2' | awk '{{total += length($0); count++}} END {{print count, total}}'"
	try:
		count, total_bases = map(int, subprocess.check_output(cmd, shell=True).strip().split())
		return count, total_bases
	except subprocess.CalledProcessError as e:
		logging.error(f"Error reading {fastq}: {e}")
		return 0, 0

def count_contigs(fasta):
	"""Count the number of contigs and their total size in a gzipped FASTA file."""
	if not fasta or not os.path.isfile(fasta):
		return 0, 0
	cmd = f"zcat {fasta} | awk '/^>/ {{if (seqlen) print seqlen; seqlen=0; next;}} {{seqlen += length($0);}} END {{if (seqlen) print seqlen;}}'"
	try:
		contig_lengths = list(map(int, subprocess.check_output(cmd, shell=True).strip().split()))
		return len(contig_lengths), sum(contig_lengths)
	except subprocess.CalledProcessError as e:
		logging.error(f"Error reading {fasta}: {e}")
		return 0, 0

def count_n_bases(fasta):
	"""Count the number of 'N' bases in a gzipped FASTA file."""
	if not fasta or not os.path.isfile(fasta):
		return 0
	cmd = f"zcat {fasta} | grep -v '^>' | tr -cd 'Nn' | wc -c"
	try:
		return int(subprocess.check_output(cmd, shell=True).strip())
	except subprocess.CalledProcessError as e:
		logging.error(f"Error reading {fasta}: {e}")
		return 0


def trim_and_assembly(row):
	"""
	Process trimming, assembly, and statistics generation for a given BioSample and Run.
	"""
	# Define paths
	sample = f"{row['BioSample']}.{row['Run']}"
	assembly = os.path.join(config['assemblies'], f"{sample}.fna.gz")
	denovo_path = os.path.join(config['denovo'], f"{sample}_denovo")
	statfile = os.path.join(config['stats'], f"{sample}.stats.tsv")
	# Step 1: Check if assembly exists and is not corrupt
	if os.path.isfile(assembly) and check_gz_integrity(assembly):
		msg(f"{assembly} is ready, skipping processing.")
	else:
		# Step 2: Download raw reads if necessary
		if not download_raw_reads(row, sample):
			return

		# Step 3: Trim raw reads if necessary
		if not trim_reads(row, sample):
			return  # Skip if trimming fails

		# Step 4: Perform assembly using SPAdes if necessary
		if not perform_assembly(row, sample, denovo_path):
			return  # Skip if assembly fails

		# Step 5: Filter contigs
		msg(f"{sample}: Filtering contigs...")
		if not collect_scaffolds(sample, assembly, denovo_path):
			return #Break if we fail at this point

		# Step 6: Generate statistics
		generate_statistics(row, sample, assembly, statfile)

	# Step 7: Remove intermediate files - only if everything was successful
	cleanup_intermediate_files(row, assembly, denovo_path, deletion=True)

	msg(f"{sample}: Processing complete.")


def download_raw_reads(row, sample):
	"""Download raw reads if they do not exist."""
	msg(f"Downloading {sample}")
	status = {item: not item or (item and os.path.isfile(item) and check_gz_integrity(item)) for item in [row["SE"], row["R1"], row["R2"]]}
	if all(status):
		msg(f"{sample}: Raw reads are already downloaded.")
		return True
	else:
		msg(f"{sample}: Downloading raw reads...")
		print(f"fasterq-dump {row['Run']} -q -t {config['fasterq']} -o {config['fasterq']} -m 1GB")
		exit()
		for item in [row['SE'], row['R1'], row['R2']]:
			if item:
				if os.path.isfile(source := os.path.join(config['fasterq'], Path(item).stem)):
					sysexec(f"gzip -f {source}")
					shutil.copy(f"{source}.gz", item)
	for item in [row["SE"], row["R1"], row["R2"]]:
		status[item] = not item or (item and os.path.isfile(item) and check_gz_integrity(item))
	if not all([status[item] for item in [row["SE"], row["R1"], row["R2"]]]):
		msg(f"{sample}: Error with downloading, missing files")
		return False
		


def trim_reads(row, sample):
	"""Trim raw reads if they are not already trimmed."""
	if not all(not item or (item and os.path.isfile(item)) for item in [row["SE trim"], row["R1 trim"], row["R2 trim"]]):
		msg(f"{sample}: Checking and trimming raw reads...")
		if not check_set_of_fq(row["Run"], row["SE"], row["R1"], row["R2"]):
			error(f"{sample}: Invalid raw fastq files, skipping...")
			return False
		msg(f"{sample}: Performing trimming...")
		sysexec(row["trim"])
	else:
		msg(f"{sample}: Reads are already trimmed.")
	return True


def perform_assembly(row, sample, denovo_path):
	"""Run SPAdes assembly if not already done."""
	scaffolds_path = os.path.join(denovo_path, "scaffolds.fasta")
	if os.path.isfile(scaffolds_path):
		msg(f"{sample}: Assembly already completed, skipping SPAdes.")
		return True

	msg(f"{sample}: Verifying trimmed reads...")
	if not check_set_of_fq(f"{row['Run']} (trimmed)", row["SE trim"], row["R1 trim"], row["R2 trim"]):
		error(f"{sample}: Invalid trimmed files, skipping...")
		return False

	msg(f"{sample}: Running SPAdes...")
	sysexec(row["spades"])
	return True


def generate_statistics(row, sample, assembly, statfile):
	"""Generate statistics for trimming and assembly."""
	msg(f"{sample}: Generating statistics...")
	stats = OrderedDict([
		("Run", row["Run"]),
		("BioSample", row["BioSample"]),
		("SE Reads Before Trimming", 0),
		("R1 Reads Before Trimming", 0),
		("R2 Reads Before Trimming", 0),
		("SE Reads After Trimming", 0),
		("R1 Reads After Trimming", 0),
		("R2 Reads After Trimming", 0),
		("Total Nucleotides Before Trimming", 0),
		("Total Nucleotides After Trimming", 0),
		("Percentage of Kept Nucleotides", 0),
		("N Count", 0),
		("Contig Number", 0),
		("Genome Size", 0),
		("Average Coverage", 0),
	])

	# Populate read counts and nucleotide statistics
	for fq, key_reads, key_bases in [
		(row["SE"], "SE Reads Before Trimming", "Total Nucleotides Before Trimming"),
		(row["R1"], "R1 Reads Before Trimming", "Total Nucleotides Before Trimming"),
		(row["R2"], "R2 Reads Before Trimming", "Total Nucleotides Before Trimming"),
		(row["SE trim"], "SE Reads After Trimming", "Total Nucleotides After Trimming"),
		(row["R1 trim"], "R1 Reads After Trimming", "Total Nucleotides After Trimming"),
		(row["R2 trim"], "R2 Reads After Trimming", "Total Nucleotides After Trimming"),
	]:
		read_count, total_bases = count_reads_and_bases(fq)
		stats[key_reads] += read_count
		stats[key_bases] += total_bases

	# Calculate percentages and other statistics
	if stats["Total Nucleotides Before Trimming"] > 0:
		stats["Percentage of Kept Nucleotides"] = round(
			100 * stats["Total Nucleotides After Trimming"] / stats["Total Nucleotides Before Trimming"], 2
		)
	if os.path.isfile(assembly):
		stats["Contig Number"], stats["Genome Size"] = count_contigs(assembly)
		stats["N Count"] = count_n_bases(assembly)
		if stats["Genome Size"] > 0:
			stats["Average Coverage"] = stats["Total Nucleotides After Trimming"] / stats["Genome Size"]

	# Write statistics to a TSV file
	with open(statfile, "w", newline="") as tsvfile:
		writer = csv.DictWriter(tsvfile, fieldnames=stats.keys(), delimiter="\t")
		writer.writeheader()
		writer.writerow(stats)
	msg(f"{sample}: Statistics written to {statfile}.")


def cleanup_intermediate_files(row, assembly, denovo_path, deletion=False):
	"""Delete intermediate files if the assembly exists and is valid."""
	deleted_files = []
	spades_log = os.path.join(config['denovo'], f"{row['BioSample']}.{row['Run']}.spades.log")
	spades_err = os.path.join(config['denovo'], f"{row['BioSample']}.{row['Run']}.spades.err")
	for item in [row["SE"], row["R1"], row["R2"], row["SE trim"], row["R1 trim"], row["R2 trim"], spades_err, spades_log]:
		if item and os.path.isfile(item):
			if deletion:
				os.remove(item)
			deleted_files.append(item)
	if os.path.isdir(denovo_path):
		if deletion:
			shutil.rmtree(denovo_path)
		deleted_files.append(denovo_path)
	if deleted_files:
		msg(f"{row['BioSample']}.{row['Run']}: Deleted intermediate files: {', '.join(deleted_files)}")


def process_row_with_progress(row, progress_counter, active_counter, total_rows):
	"""
	Process a single row and update progress counters.

	:param row: The data row to process.
	:param progress_counter: Shared counter for completed samples.
	:param active_counter: Shared counter for currently active samples.
	:param total_rows: Total number of rows to process.
	"""
	# Increment the active counter
	active_counter.value += 1

	# Print progress
	print(
		f"Processed: {progress_counter.value}, "
		f"Processing: {active_counter.value}, "
		f"Remaining: {total_rows - progress_counter.value - active_counter.value}",
		end="\r"
	)

	# Perform the processing
	trim_and_assembly(row)

	# Update counters after processing
	active_counter.value -= 1
	progress_counter.value += 1

	# Print updated progress
	print(
		f"Processed: {progress_counter.value}, "
		f"Processing: {active_counter.value}, "
		f"Remaining: {total_rows - progress_counter.value - active_counter.value}",
		end="\r"
	)



def run_trim_and_assembly_for_all_samples(commands):
	"""
	Run trimming and assembly for all samples using multiprocessing.

	:param commands: Path to the commands TSV file.
	"""
	# Step 1: Create required directories
	for item in ['trimmed', 'denovo', 'assemblies', 'stats', 'fasterq']:
		mkdir_force(config[item])
		
	assembled_genomes = [os.path.basename(x) for x in glob.glob(os.path.join(config['assemblies'],"*.fna.gz"))]
	
	# Step 2: Load commands from the TSV file
	with open(commands, "r") as f:
		reader = csv.DictReader(f, delimiter="\t")
		rows = [row for row in reader if f"{row['BioSample']}.{row['Run']}.fna.gz" not in assembled_genomes]

	total_rows = len(rows)
	msg(f"Processing {total_rows} genomes", onscreen=True)
	
	##for row in rows:
	##	trim_and_assembly(row)
	##	exit()
	# Step 3: Set up multiprocessing.Manager to manage shared state
	with mp.Manager() as manager:
		progress_counter = manager.Value('i', 0)  # Shared counter for completed samples
		active_counter = manager.Value('i', 0)    # Shared counter for currently active samples

		# Step 4: Run trim_and_assembly in parallel
		with mp.Pool(pargs.n) as pool:  # pargs.n specifies the number of processes
			pool.starmap(
				process_row_with_progress,
				[(row, progress_counter, active_counter, total_rows) for row in rows]
			)

	# Step 5: Final message
	info("Assembly completed for all samples.")

	
	
def main():
	"""Main execution flow"""
	global config
	global pargs
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
	####download_reads(config['urinary_reads'], config["urinary_with_dates_and_locations"], config["raw_reads"])
	#
	get_reference_genome()
	#
	commands = prepare_assembly_commands(config["urinary_with_dates_and_locations"], config['urinary_reads'], config['commands'])
	#This is the part for trimming and assembly
	run_trim_and_assembly_for_all_samples(config['commands'])




	
if __name__ == "__main__":
	init_logging()
	main()
