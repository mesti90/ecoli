`python3 ecoli/bin/download_species_from_ncbi.py3`

Downloads all genomes found in NCBI refseq and genbank belonging to the specified species (E. coli)

Also downloads related biosample metadata

`python3 ecoli/bin/fetch_biosample_data.py3`

Converts downloaded biosample metadata to table

`python3 ecoli/bin/extract_urinary_genomes.py3`

Selects biosamples which have geolocation and collection date, human origin, and which match any of the keywords:

"CAUTI", "periurethral", "uroepithelium", "uroepithelial", "renal", "tubular epithelial", "umbrella cells", "superficial facet cells", "uroplakin", "uropathogenic", "ureter", "urination", "prostatitis", "pyelonephritis", "pyrexia", "Genitourinary", "Urinary", "Urine", "UTI", "Bladder", "Cystitis", "Catheter", "UPEC", "kidney", "urethra", "urethral", "nephron", "Nephritis"
