#!/usr/bin/python3
'''Automatically downloading all genomes belonging to a certain species from NCBI refseq and genbank databases
'''

import csv
import os
import argparse
import multiprocessing as mp
import time


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-dd",help="Folder of NCBI metadata",default="NCBI_data")
parser.add_argument("-s",help="Species",default="Escherichia coli")
parser.add_argument("-n",help="Number of threads",default=80,type=int)
parser.add_argument("-gd",help="Directory of genomes",default="Genomes")
pargs = parser.parse_args()


NCBI_DIR = pargs.dd
GENOME_DIR = pargs.gd
REFSEQ = os.path.join(NCBI_DIR,"assembly_summary_refseq.txt")
GENBANK = os.path.join(NCBI_DIR,"assembly_summary_genbank.txt")
SPECIES = "_".join(pargs.s.split())
SELECTED_REFSEQ = os.path.join(NCBI_DIR,f"assembly_summary_refseq_{SPECIES}.txt")
SELECTED_GENBANK = os.path.join(NCBI_DIR,f"assembly_summary_genbank_{SPECIES}.txt")
GENOMES = os.path.join(NCBI_DIR,f"{SPECIES}.txt")
DOWNLOAD = os.path.join(NCBI_DIR,"download_commands.txt")
GENOMES_DOWNLOADED = os.path.join(NCBI_DIR,"genomes_downloaded.flag")
UNIQUE_GENOMES = os.path.join(NCBI_DIR,f"{SPECIES}.unique_biosample.tsv")
BIOSAMPLE_DATA = os.path.join(NCBI_DIR,f"{SPECIES}.biosample.data.txt")
FAILED = os.path.join(NCBI_DIR,f"{SPECIES}.biosample.failed.log")
EFETCH_LOG = os.path.join(NCBI_DIR,f"{SPECIES}.biosample.download.log")

def mkdir_force(path):
	os.makedirs(path, exist_ok = True)

def sysexec(cmd):
	msg(cmd)
	os.system(cmd)

def wget(link):
	sysexec(f"wget {link}")


def download_assembly_summaries():
	if not os.path.isfile(REFSEQ):
		wget(f"-nc https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O {REFSEQ}")
	else:
		print(f"{REFSEQ} exists, I refuse to overwrite it. If you want to download it again, please remove it on your own.")
	if not os.path.isfile(GENBANK):
		wget(f"-nc https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt -O {GENBANK}")
	else:
		print(f"{REFSEQ} exists, I refuse to overwrite it. If you want to download it again, please remove it on your own.")
	msg(f"Creating {GENOMES}")
	if os.path.isfile(GENOMES):
		msg("f{GENOMES} is already created")
		return
		#os.remove(GENOMES)
	filter_summary(REFSEQ,GENOMES)
	filter_summary(GENBANK,GENOMES)
	msg(f"Ready with {GENOMES}")


def filter_summary(source,result):
	with open(source) as f, open(result,"a") as g:
		rdr = csv.reader(f,delimiter='\t')
		wtr = csv.writer(g,delimiter="\t")
		for row in rdr:
			if len(row) > 7:
				block = row[7].strip("[").strip("]").split()
				if " ".join(block[:2]) == pargs.s:
					wtr.writerow(row)
	#sysexec(fr'''awk -F '\t' '{{if (($8 == "Klebsiella pneumoniae" )) {{print}}}}' {source} >> {result}''')

def smartint(x):
	try:
		return int(x)
	except ValueError:
		return 0
	except:
		return -1

def collect_unique_identifiers(genomes,download_file,unique_genomes_file):
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

def msg(cmd):
	uzenet = f'{time.strftime("%H:%M:%S")} {cmd}'
	print(uzenet)



def read_file_lines(fname,split=False):
	msg(f"Reading {fname}")
	if not os.path.isfile(fname):
		return []
	lines = []
	with open(fname) as f:
		for line in f:
			line = line.strip()
			if len(line) and line[0] != "#":
				if split:
					lines.append(line.split()[0])
				else:
					lines.append(line)
	return list(set(lines))


def download_all(dl_file,downloaded):
	if os.path.isfile(downloaded):
		print("Genomes are already downloaded, skipping")
		return
	pool = mp.Pool(pargs.n)
	links = read_file_lines(dl_file)
	pool.map(download,links)
	sysexec(f"touch {downloaded}")

def download(link):
	result = link.split()[-1].rsplit(".",1)[0]
	if os.path.isfile(result):
		return
	sysexec(link)

def extract_biosample_data(genomes,output,failed,log):
	sysexec(f'''tail -n +2 {genomes} | cut -f 2 | xargs -P 8 -I {{}} sh -c 'efetch -db biosample -format full -id {{}} || echo "Failed ID: {{}}" >> {failed}' > {output} 2> {log}''')

def main():
	mkdir_force(NCBI_DIR)
	mkdir_force(GENOME_DIR)
	download_assembly_summaries()
	collect_unique_identifiers(GENOMES,DOWNLOAD,UNIQUE_GENOMES)
	download_all(DOWNLOAD,GENOMES_DOWNLOADED)
	extract_biosample_data(UNIQUE_GENOMES,BIOSAMPLE_DATA,FAILED,EFETCH_LOG)


if __name__ == "__main__":
	main()
