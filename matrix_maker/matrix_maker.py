#! /usr/bin/python
"""
matrix_maker.py

A simple script that finds the NCBI taxid for a list of taxa, and
then searches GenBank for a given gene using the taxids. The script
then downloads the sequences and produces a MAFFT alignment.

Copyright 2015 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


# Configurations: edit these!
email = "freyman@berkeley.edu"
gene_names = ["rbcL", "rbcl", "RBCL", "Rbcl"]
taxa_file = "input.txt"


def get_sequences(taxid):
    """
    Searches Entrez Nucleotide database for taxid and gene names and 
    downloads results.
    """
    from Bio import Entrez
    Entrez.email = email
    # (txid202994[Organism] AND (rbcL[All Fields] OR internal transcribed spacer[All Fields])
    term = "txid" + taxid + "[Organism] AND ("
    for i, gene in enumerate(gene_names):
        if i == 0:
            term = term + gene + "[All Fields]"
        else:
            term = term + " OR " + gene + "[All Fields]"
    term = term + ")"
    print("Using search term: " + term)
    handle = Entrez.esearch(db="nuccore", term=term)
    records = Entrez.read(handle)
    gi_list = records["IdList"]
    gi_str = ",".join(gi_list)
    print("Found GenBank GIs: " + gi_str)
    handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
    from Bio import SeqIO
    records = SeqIO.parse(handle, "gb")
    final_records = []
    for record in records:
        final_records.append(record)
    return final_records


def get_taxon_id(taxon):
    """
    Gets the taxid from entrez taxonomy.
    """
    import urllib
    import re
    toolname = "matrix_maker"
    params = {
        'db': 'taxonomy',
        'tool': toolname,
        'email': email,
        'term': taxon,
        'rettype': 'xml',
    }
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
    url = url + urllib.urlencode(params)
    data = urllib.urlopen(url).read()
    if re.search('<Id>(\S+)</Id>', data):
        taxid = re.search('<Id>(\S+)</Id>', data).group(1)
    else:
        taxid='not found'
    return taxid


def main():
    
    print("\n\nmatrix_maker.py\n\n")
    
    print("Getting all taxid...\n")
    name_file = open(taxa_file)
    names = name_file.readlines()
    taxids = []
    import time
    for name in names:
        name = "%s" %(name.split()[0])
        taxid = get_taxon_id(name)
        print(name + "\t" + taxid)
        taxids.append( taxid )
        # dont overload genbank
        time.sleep(0.1)

    print("\nDownloading sequences for each taxid...\nKeeping the longest sequence for each taxon...\n")
    from Bio import Entrez
    from Bio import SeqIO
    final_records = []
    for taxid in taxids:
        if taxid != "not found":
            records = get_sequences(taxid)
            # find the longest sequence
            longest_len = 0
            longest_seq = None
            for record in records:
                if len(record) > longest_len:
                    longest_len = len(record)
                    longest_seq = record
            if longest_seq != None:
                final_records.append(longest_seq)
    
    print("\nGenerating unaligned FASTA file...\n")
    SeqIO.write(final_records, "output_unaligned.fasta", "fasta")

    print("Making alignment with MAFFT...")
    from Bio.Align.Applications import MafftCommandline
    mafft_cline = MafftCommandline(input="output_unaligned.fasta")
    mafft_cline.set_parameter("--auto", True)
    mafft_cline.set_parameter("--adjustdirection", True)
    print(str(mafft_cline))
    stdout, stderr = mafft_cline()
    print("Writing alignment to FASTA file...\n")
    with open("output_aligned.fasta", "w") as handle:
        handle.write(stdout)
    print("Done!\n")


if __name__ == "__main__":
    main()
