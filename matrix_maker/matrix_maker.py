#! /usr/bin/python
"""
matrix_maker.py

A simple script that finds the NCBI taxid for a list of taxa, and
then searches GenBank for a given gene using the taxids. The script
then downloads the sequences and produces a MAFFT alignment.

Copyright 2015 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


import argparse
import time
from taxon import Taxon


# Configurations: edit these!
#email = "freyman@berkeley.edu"
#gene_names = ["rbcL", "rbcl", "RBCL", "Rbcl"]
#taxa_file = "input.txt"



def main():

    # parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--email", "-e", help="Email address for NCBI database searches.")
    parser.add_argument("--genes", "-g", help="Text file that contains a list of all gene names.")
    parser.add_argument("--species", "-s", help="Text file that contains a list of all species binomials and their synonyms.")
    parser.add_argument("--taxids", "-t", help="Text file that contains a list of all taxids. Optional, use this to avoid repeating the NCBI taxid lookups.")
    args = parser.parse_args()

    print("\n\nmatrix_maker.py\n\n")

    if not args.email:
        print("NCBI requires an email address for database searches. Please use the --email flag to specify an email address.")
        sys.exit(0)

    if not args.species:
        print("Please specify a list of taxa to search for.")
        sys.exit(0)

    if not args.genes:
        print("Please specify a list of genes to search for.")
        sys.exit(0)

    # list of all taxon objects
    taxa = []

    # check for taxid



    print("Getting all taxid...\n")
    print("Writing taxids to file taxids.txt...\n")
    taxids_file = open("taxids.txt", "w")
    name_file = open(taxa_file)
    names = name_file.readlines()
    taxids = []
    for name in names:
        name = "%s" %(name.split()[0])
        taxid = get_taxon_id(name)
        name_taxid_text = name + "\t" + taxid
        print(name_taxid_text)
        taxids_file.write(name_taxid_text + "\n")
        taxids.append( taxid )
        # dont overload genbank
        time.sleep(0.1)
    taxids_file.close()

    print("\nDownloading sequences for each taxid...\n") #Keeping the longest sequence for each taxon...\n")
    from Bio import Entrez
    from Bio import SeqIO
    final_records = []
    for taxid in taxids:
        if taxid != "not found":
            records = get_sequences(taxid)
            # keep all records
            final_records = final_records + records
            # dont overload genbank
            time.sleep(0.2)

            # find the longest sequence
            #longest_len = 0
            #longest_seq = None
            #for record in records:
            #    if len(record) > longest_len:
            #        longest_len = len(record)
            #        longest_seq = record
            #if longest_seq != None:
            #    final_records.append(longest_seq)
    
    print("\nGenerating unaligned FASTA file with GenBank formatted description...\n")
    SeqIO.write(final_records, "output_unaligned_gb_format.fasta", "fasta")

    print("Generating unaligned FASTA file with custom formatted description...\n")
    unaligned_file = open("output_unaligned_custom_format.fasta", "w")
    for record in final_records:
        # remove the organism name from the description
        description = record.description
        if description.find(record.annotations["organism"] + " ") != -1:
            description = description.replace(record.annotations["organism"] + " ", "")
        # custom format for Andrew: >Organism name_accession_description
        description = record.annotations["organism"] + "_" + record.id + "_" + description
        description = description.replace(" ", "_")
        unaligned_file.write(">" + description + "\n")
        unaligned_file.write(str(record.seq) + "\n")
    unaligned_file.close()


#    print("Making alignment with MAFFT...")
#    try:
#        from Bio.Align.Applications import MafftCommandline
#        mafft_cline = MafftCommandline(input="output_unaligned_custom_format.fasta")
#        mafft_cline.set_parameter("--auto", True)
#        mafft_cline.set_parameter("--adjustdirection", True)
#        print(str(mafft_cline))
#        stdout, stderr = mafft_cline()
#        print("Writing alignment to FASTA file...\n")
#        with open("output_aligned.fasta", "w") as handle:
#            handle.write(stdout)
#    except:
#        print("Problem finding MAFFT, alignment skipped.")
        
    print("Done!\n")


if __name__ == "__main__":
    main()
