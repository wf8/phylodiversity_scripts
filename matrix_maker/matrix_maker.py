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
import csv
import os
import sys
import time
from Bio import SeqIO
from taxon import Gene
from taxon import Taxon


def main():

    # parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--email", "-e", help="Email address for NCBI database searches.")
    parser.add_argument("--genes", "-g", help="Text file that contains a list of all gene names.")
    parser.add_argument("--max_seq_length", "-m", help="Optional. Sets the maximum sequence length to include. Use this to exclude genomes.")
    parser.add_argument("--species", "-s", help="Text file that contains a list of all species binomials and their synonyms.")
    parser.add_argument("--taxids", "-t", help="Optional. Text file that contains a list of all taxids. Use this to avoid repeating the NCBI taxid lookups.")
    args = parser.parse_args()

    print("\n\nmatrix_maker.py\n\n")

    if not args.email:
        print("NCBI requires an email address for database searches. Please use the --email flag to specify an email address.\n")
        sys.exit(0)
    else:
        email = args.email

    if not args.species or not os.path.isfile(args.species):
        print("Please specify a valid list of taxa to search for.\n")
        sys.exit(0)

    if args.max_seq_length:
        max_seq_length = int(args.max_seq_length)
    else:
        max_seq_length = -1

    genes = []
    if not args.genes or not os.path.isfile(args.genes):
        print("Please specify a valid list of genes to search for.\n")
        sys.exit(0)
    else:
        # read in gene names....
        # format of file:
        # gene_name,include,rbcL,RBCL
        # gene_name,exclude,RRRBCL
        with open(args.genes, 'rb') as csvfile:
            genereader = csv.reader(csvfile, delimiter=",")
            for row in genereader:
                if row[1] == "include":
                    gene = Gene(row[0])
                    for i in range(2, len(row)):
                        if row[i] != "":
                            gene.gene_names.append(row[i])
                    genes.append(gene)
                if row[1] == "exclude":
                    for gene in genes:
                        if gene.name == row[0]:
                            for i in range(2, len(row)):
                                if row[i] != "":
                                    gene.exclusions.append(row[i])
                    
    # list of all taxon objects
    taxa = []

    # check for taxid
    print("Checking for taxids text file...")
    if args.taxids and os.path.isfile(args.taxids):
        with open(args.taxids, 'rb') as csvfile:
            print("Found taxids text file, reading taxids...\n")
            taxidsreader = csv.reader(csvfile, delimiter=",")
            for row in taxidsreader:
                taxa.append(taxon(row[0], row[1]))
    else:
        print("No taxids text file found.\n")

    # open species list file, get synonyms and any missing taxids
    with open(args.species, 'rb') as csvfile:
        print("Checking list of species, getting missing taxids from NCBI...")
        taxids_file = open("taxids.txt", "w")
        namesreader = csv.reader(csvfile, delimiter=",")
        i = 1
        num_lines = sum(1 for line in open(args.species))
        for row in namesreader:
            # update status
            percent = str(round(100 * i/float(num_lines), 2))
            sys.stdout.write('\r' + 'Completed: ' + str(i) + '/' + str(num_lines) + ' (' + percent + '%)')
            i += 1
            # check to see if we already have a taxid for this species
            found = False
            for taxon in taxa:
                if taxon.binomial == row[0]:
                    found = True
                    taxids_file.write(taxon.binomial + "," + taxon.taxid + "\n")
                    # add synonyms
                    for i in range(1, len(row)):
                        taxon.synonyms.append(row[i])
                    break
            if not found:
                # get the taxid from NCBI
                taxon = Taxon(row[0])
                taxon.get_taxid(email)
                # dont overload genbank
                time.sleep(0.1)
                taxids_file.write(taxon.binomial + "," + taxon.taxid + "\n")
                # add synonyms
                for i in range(1, len(row)):
                    taxon.synonyms.append(row[i])
                taxa.append(taxon)
        taxids_file.close()
        print("\nWriting all taxids to file taxids.txt...")

    print("\nDownloading sequences from NCBI...\n") 
    for gene in genes:
        print("Searching for gene: " + gene.name)
        i = 1
        for taxon in taxa:
            # update status
            percent = str(round(100 * i/float(len(taxa)), 2))
            sys.stdout.write('\r' + 'Completed: ' + str(i) + '/' + str(num_lines) + ' (' + percent + '%)')
            i += 1
            if taxon.taxid != "not found":
                taxon.get_sequences(email, gene)
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

        print("\nGenerating unaligned FASTA file...\n")
        unaligned_file = open(gene.name + ".fasta", "w")
        for taxon in taxa:
            #TODO: optionally include only the longest sequence
            for record in taxon.sequences:
                if len(record) < max_seq_length and max_seq_length != -1: 
                    # custom format for Andrew: >Organism name_accession_description
                    description = taxon.binomial + "_" + record.id + "_" + record.description
                    description = description.replace(" ", "_")
                    unaligned_file.write(">" + description + "\n")
                    unaligned_file.write(str(record.seq) + "\n")
        unaligned_file.close()
    print("Done!\n")


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
        


if __name__ == "__main__":
    main()
