"""
Copyright 2015 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


import urllib
import re
from Bio import Entrez
from Bio import SeqIO



class Gene(object):
    """
    Class responsible for managing the search terms of genes we are searching for in NCBI.
    """


    name = ""           # name of gene region
    gene_names = []     # list of search terms for this gene
    exclusions = []     # list of search terms to exclude for this gene

    
    def __init__(self, name=""):
        self.name = name



class Taxon(object):
    """
    Class responsible for managing the data for each taxon.
    """


    binomial = ""           # genus_species
    taxid = ""              # NCBI taxid
    synonyms = []           # list of synonym
    sequences = []          # list of lists of Bio.SeqRecord (each list is a different gene region)


    def __init__(self, binomial, taxid=""):
        """
        Optionally accept the NCBI taxid.
        """
        self.binomial = binomial
        self.taxid = taxid


    def get_taxid(self, email):
        """
        Gets the NCBI taxid from entrez taxonomy.
        """
        toolname = "matrix_maker"
        params = {
            'db': 'taxonomy',
            'tool': toolname,
            'email': email,
            'term': self.binomial,
            'rettype': 'xml',
        }
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
        url = url + urllib.urlencode(params)
        data = urllib.urlopen(url).read()
        if re.search('<Id>(\S+)</Id>', data):
            self.taxid = re.search('<Id>(\S+)</Id>', data).group(1)
        else:
            # taxid was not found using the binomial,
            # so now check for synonyms...
            for synonym in self.synonyms:
                params = {
                    'db': 'taxonomy',
                    'tool': toolname,
                    'email': email,
                    'term': synonym,
                    'rettype': 'xml',
                }
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
                url = url + urllib.urlencode(params)
                data = urllib.urlopen(url).read()
                if re.search('<Id>(\S+)</Id>', data):
                    self.taxid = re.search('<Id>(\S+)</Id>', data).group(1)
                    break
            # if taxid is still not found
            if self.taxid == '':
                self.taxid='not found'
        return self.taxid


    def get_sequences(self, email, gene):
        """
        Searches Entrez Nucleotide database for taxid and and a list of gene names and
        downloads results. Appends the resulting list of Bio.SeqRecords to self.sequences.
        """
        Entrez.email = email
        # TODO: handle search term exclusions!!!
        # (txid202994[Organism] AND (rbcL[All Fields] OR internal transcribed spacer[All Fields])
        term = "txid" + self.taxid + "[Organism] AND ("
        for i, name in enumerate(gene.gene_names):
            if i == 0:
                term = term + name + "[All Fields]"
            else:
                term = term + " OR " + name + "[All Fields]"
        term = term + ")"
        #print("Using search term: " + term)
        handle = Entrez.esearch(db="nuccore", term=term)
        records = Entrez.read(handle)
        gi_list = records["IdList"]
        gi_str = ",".join(gi_list)
        #print("Found GenBank GIs: " + gi_str)
        handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
        records = SeqIO.parse(handle, "gb")
        self.sequences.append(records)
        return self.sequences
