#! /usr/bin/python
"""
Replace all occurences in input file of first column of csv file 
with second column and write output file.
"""

input_file = "Eucalypt_Final_PD.tre"
output_file = "Eucalypt_Final_PD_new_tips.tre"
csv_file = "Test_taxonomy_2.csv" 

import csv

with open(input_file, "r") as f:
    s = f.read()
    with open(csv_file, "rb") as f:
        csvreader = csv.reader(f)
        for row in csvreader:
            s = s.replace(row[0], row[1])
        with open(output_file, "w") as o:
            o.write(s)


