#!/usr/bin/python3

import csv
import pandas as pd

from collections import defaultdict
from functools import reduce


fname = "biosample_data.csv"
selection = "biosample_data.selected_columns.csv"
##related_biosamples = "biosamples.urinary.csv"
related_biosamples = "biosamples.urinary.csv"
related_biosamples = "biosamples.urinary.csv"
related_biosamples_nonempty = "biosamples.urinary.nonempty.csv"

def count_columns_and_rows(fname):
	print(f"Reading {fname}")
	with open(fname) as f:
		rdr = csv.reader(f, delimiter='\t')
		header = next(rdr)
		counts = {colname:0 for colname in header}
		rownum = 0
		for row in rdr:
			rownum += 1
			for item,colname in zip(row,header):
				if item:
					counts[colname] += 1
	header.sort(key=lambda x:-counts[x])
	outf = f"{fname}.counts"
	print(f"Writing {outf}")
	with open(outf,"w") as g:
		wtr = csv.writer(g,delimiter='\t')
		wtr.writerow(["Attribute","Total","Percentage"])
		wtr.writerow(["Sample number",rownum,100.00])
		for colname in header:
			percentage = counts[colname] / rownum * 100 if rownum > 0 else 0
			wtr.writerow([colname, counts[colname], f"{percentage:.02f}"])



#We keep only those samples where "geographic location" and "collection date" are non empty
def select_columns(input_file, output_file, columns_to_select):
	with open(input_file, newline='') as infile, open(output_file, mode='w', newline='') as outfile:
		rdr = csv.DictReader(infile, delimiter="\t")
		# Filter for columns that exist in the input file
		selected_columns = [col for col in columns_to_select if col in rdr.fieldnames]
		writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=selected_columns)
		writer.writeheader()
		for row in rdr:
			selected_row = {col: row[col] for col in selected_columns}
			if row['geographic location'] and row['collection date'] and row['collection date'] != "missing" and row['geographic location'] != "missing" and row['Organism'].startswith("Escherichia coli"):
				writer.writerow(selected_row)

def filter_rows(input_file, output_file, keywords, human_keywords):
	with open(input_file) as f, open(output_file, 'w') as g:
		g.write(next(f))
		for line in f:
			lower_line = line.lower()
			if any(keyword.lower() in lower_line for keyword in keywords):#! and any(kw.lower() in lower_line for kw in human_keywords) :
				g.write(line)


def drop_empty_columns_and_unrelated_genomes(input_file, output_file):
	df = pd.read_csv(input_file, sep="\t", low_memory=False)
	
	conditions = [
	df['Organism'].str.startswith("Escherichia coli"),
	df['collection date'].str.strip().ne(""),
	df['collection date'].str.strip().ne("missing"),
	df['geographic location'].str.strip().ne(""),
	df['geographic location'].str.strip().ne("missing")
	]
	combined_condition = reduce(lambda x, y: x & y, conditions)


	
	df_filtered = df[combined_condition].copy()
	df_cleaned = df_filtered.dropna(axis=1, how='all')
	df_cleaned.to_csv(output_file, index=False, sep="\t")

columns_to_select = ["BioSample", "collection date", "geographic location",  "isolation source", "strain", "source type", "host", "host description", "host disease", "host disease outcome", "study disease", "serotype", "pathotype", "environmental sample", "sample type", "Entry notes"]
	
keywords = ["CAUTI", "periurethral", "uroepithelium", "uroepithelial", "renal", "tubular epithelial", "umbrella cells", "superficial facet cells", "uroplakin", "uropathogenic", "ureter", "urination", "prostatitis", "pyelonephritis", "pyrexia", "Genitourinary", "Urinary", "Urine", "UTI", "Bladder", "Cystitis", "Catheter", "UPEC", "kidney", "urethra", "urethral", "nephron", "Nephritis"]

human_keywords = ["Homo", "homo", "human"]





if __name__ == "__main__":
	os.chdir("NCBI_data")
	count_columns_and_rows(fname)
	####!select_columns(fname, selection, columns_to_select)
	####!filter_rows(selection, related_biosamples, keywords, human_keywords)
	filter_rows(fname, related_biosamples, keywords, human_keywords)
	drop_empty_columns_and_unrelated_genomes(related_biosamples, related_biosamples_nonempty)
	os.chdir("..")
