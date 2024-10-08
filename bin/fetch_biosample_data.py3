#!/usr/bin/python3
import re
import csv
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", default="Escherichia_coli.biosample.data.txt", help="Path with metadata retrieved from NCBI")
parser.add_argument("-o", default="biosample_data.csv", help="Output file")
pargs = parser.parse_args()


# Parse raw text and extract records
def parse_ncbi_metadata(fname):
	entry_num = 0
	entries = []
	my_entry = defaultdict(str)
	with open(fname) as f:
		for line in f:
			line = line.strip()
			if line.startswith("1: "):
				entry_num += 1
				if my_entry:
					entries.append(my_entry)
				my_entry = defaultdict(str)
				my_entry["Entry notes"] = line[3:]
			elif line.startswith("Organism: "):
				block = line.split(": ")
				my_entry['Organism'] = block[1]
			elif line.startswith("Identifiers: "):
				id_list = line.replace("Identifiers: ","").split(";")
				for item in id_list:
					block = item.strip().split(": ")
					if len(block) == 2:
						my_entry[block[0]] = block[1]
			elif line.startswith("/"):
				key,value = line.split("=",1)
				my_entry[key.strip("/")] = value.strip('"')
	return entries

# Write the parsed data into a CSV file
def write_to_csv(records, output_file):
	# Get all unique attribute keys
	all_attributes = set()
	for record in records:
		all_attributes.update(record.keys())
	
	first_cols = ['BioSample', 'Sample name', 'Organism']
	
	all_attributes -= set(first_cols)
	
	# Define column headers
	header = first_cols + sorted(all_attributes)
	
	# Write to CSV
	with open(output_file, 'w', newline='') as csvfile:
		writer = csv.DictWriter(csvfile, fieldnames=header, delimiter="\t")
		writer.writeheader()
		for record in records:
			row = {attribute: record.get(attribute,'') for attribute in header}
			writer.writerow(row)

if __name__ == "__main__":
	os.chdir("NCBI_data")
	records = parse_ncbi_metadata(pargs.i)
	#Write the results to CSV
	write_to_csv(records, pargs.o)
	os.chdir("..")

print(f"Data successfully written to {pargs.o}")
