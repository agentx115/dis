# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 18:04:35 2019

@author: Rhiannon
"""

#%%
import os
os.getcwd()
os.chdir("C:\\Users\\Rhiannon\\Documents\\Year_4_Uni\\project\\python")

#%%
from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
#%% Open excel file
import csv
#%%
#open hgt-lists
rmag_genes = SeqIO.parse(open("R.mag.GDRE01.1.fsa_nt"), "fasta")
rsoc_genes = SeqIO.parse(open("R.soc.GDRD01.1.fsa_nt"), "fasta")

#%%
#make file for new ones in fasta
rmag_hgt_list = []
rsoc_hgt_list = []

#%%
with open("geneID.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        #print(row["query.id"])
        rowstr = str(row["query.id"])
        if "RMAG" in rowstr:
            rmag_hgt_list.append(rowstr)
        elif "RSOC" in rowstr:
            rsoc_hgt_list.append(rowstr)
            
#%%
protein_file = open("rmag_extract.txt","w+")  
          
for fasta_rmag in rmag_genes:
    name, sequence = fasta_rmag.description, str(fasta_rmag.seq)
    #print(name)
    for gene_name in rmag_hgt_list:
        if str(gene_name) in str(name):
            print("found ", str(gene_name), " in ", str(name))
            fasta_seq = str(">" + str(name) + "\n" + str(sequence) + "\n")
            protein_file.write(fasta_seq)


protein_file.close()

#%%
protein_file_rsoc = open("rsoc_extract.txt","w+")  
          
for fasta_rsoc in rsoc_genes:
    name, sequence = fasta_rsoc.description, str(fasta_rsoc.seq)
    #print(name)
    for gene_name in rsoc_hgt_list:
        if str(gene_name) in str(name):
            print("found ", str(gene_name), " in ", str(name))
            fasta_seq = str(">" + str(name) + "\n" + str(sequence) + "\n")
            protein_file_rsoc.write(fasta_seq)


protein_file_rsoc.close()
            
            
            
            
            
            
            
            