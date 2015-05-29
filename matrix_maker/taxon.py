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


    def __init__(self, name=""):
        self.name = name         # name of gene region
        self.gene_names = []     # list of search terms for this gene
        self.exclusions = []     # list of search terms to exclude for this gene



class Taxon(object):
    """
    Class responsible for managing the data for each taxon.
    """


    def __init__(self, binomial, taxid=""):
        """
        Optionally accept the NCBI taxid.
        """
        self.binomial = binomial    # genus_species
        self.taxid = taxid          # NCBI taxid
        self.synonyms = []          # list of synonym
        self.sequences = {}         # a dictionary holding lists of Bio.SeqRecords, each key is a gene name


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
        Exclude any search term exclusions found in the records description.
        """
        Entrez.email = email
        # (txid202994[Organism] AND (rbcL[All Fields] OR internal transcribed spacer[All Fields])
        term = "txid" + self.taxid + "[Organism] AND ("
        for i, name in enumerate(gene.gene_names):
            if i == 0:
                term = term + name + "[All Fields]"
            else:
                term = term + " OR " + name + "[All Fields]"
        term = term + ")"
        #print("Using search term: " + term)
        try:
            handle = Entrez.esearch(db="nuccore", term=term)
            records = Entrez.read(handle)
            gi_list = records["IdList"]
            gi_str = ",".join(gi_list)
            #print("Found GenBank GIs: " + gi_str)
            handle = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
            records = SeqIO.parse(handle, "gb")
            final_records = []
            for record in records:
                exclude = False
                # check to make sure exclusions are not in the description
                for exclusion in gene.exclusions:
                    if exclusion in record.description:
                        exclude = True
                        break
                # check to make sure the search terms are actually in the description 
                # (they might have been in other parts of the record)
                include = False
                for name in gene.gene_names:
                    if name in record.description:
                        include = True
                        break
                if not exclude and include:
                    final_records.append(record)
        except Exception as e:
            print("\nCaught error while downloading from NCBI: " + str(e))
            print("Trying again...")
            time.sleep(0.5)
            return self.get_sequences(email, gene)
        self.sequences[gene.name] = final_records
        return final_records


    def get_longest_seq(self, gene_name, max_seq_length = -1):
        """
        Returns the longest SeqRecord shorter than max_seq_length for a given gene.
        """
        try:
            records = self.sequences[gene_name]
        except:
            records = []
        longest_len = 0
        longest_seq = None
        for record in records:
            if len(record) < max_seq_length or max_seq_length == -1:
                if len(record) > longest_len:
                    longest_len = len(record)
                    longest_seq = record
        return longest_seq
