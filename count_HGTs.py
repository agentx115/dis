# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 17:01:34 2019

@author: Rhiannon
"""

#%%
import csv
import os
os.getcwd()
os.chdir("C:\\Users\\Rhiannon\\Documents\\Year_4_Uni\\project\\python")

#%% get the separate rows

from Bio import SeqIO

#open hgt-lists
rmag_genes = SeqIO.parse(open("R.mag.GDRE01.1.fsa_nt"), "fasta")
rsoc_genes = SeqIO.parse(open("R.soc.GDRD01.1.fsa_nt"), "fasta")

#make file for new ones in fasta
rmag_hgt_list = []
rsoc_hgt_list = []

with open("geneID.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        #print(row["query.id"])
        rowstr = str(row["query.id"])
        if "RMAG" in rowstr:
            rmag_hgt_list.append(rowstr)
        elif "RSOC" in rowstr:
            rsoc_hgt_list.append(rowstr)

#%% get the correct ID's
rsoc_HGT_right = []

with open("rsoc_HGT.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        rsoc_HGT_right.append(row["Name"])
            
print(len(rsoc_HGT_right))

rsoc_genes_file = SeqIO.parse(open("rsoc_extract.txt"), "fasta")

list_of_good = []
list_of_bad = []
for fasta in rsoc_genes_file:
    name, sequence = fasta.description, str(fasta.seq)
    print(name)
    split_name = name.split(" ")
    if split_name[0] in rsoc_HGT_right:
        print(split_name[0], " is ", split_name[4], " found in the list of right ones")
        list_of_good.append(split_name[4])
    else:
        print(split_name[0], " is ", split_name[4], " which is not in the right ones")
        list_of_bad.append(split_name[4])
#%%
print(len(list_of_good), len(list_of_bad))

#%% for R. magnacalcarata
rmag_HGT_right = []

with open("rmag_HGT.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        rmag_HGT_right.append(row["Name"])
            
print(len(rmag_HGT_right))

rmag_genes_file = SeqIO.parse(open("rmag_extract.txt"), "fasta")

list_of_good_rmag = []
list_of_bad_rmag = []
for fasta in rmag_genes_file:
    name, sequence = fasta.description, str(fasta.seq)
    print(name)
    split_name = name.split(" ")
    if split_name[0] in rmag_HGT_right:
        print(split_name[0], " is ", split_name[4], " found in the list of right ones")
        list_of_good_rmag.append(split_name[4])
    else:
        print(split_name[0], " is ", split_name[4], " which is not in the right ones")
        list_of_bad_rmag.append(split_name[4])
#%%
print(len(list_of_good_rmag), len(list_of_bad_rmag))

#%% look through csv for the query.id's
import pandas as pd
#%%

data = pd.read_csv("query_id_info.csv", index_col = "query.id")

GO_ann = pd.read_csv("GO_info.csv", index_col = "GO.ID")

#%%
rmag_missing = data.loc[list_of_bad_rmag]

for index, row in rmag_missing.iterrows():
    print(index, row["go.numbers"], row["annotation"], row["best.HGT.taxon"], 
          row["second.best.HGT.taxon"], "\n")
    
rsoc_data = rmag_missing.filter(["go.numbers", "annotation", "best.HGT.taxon", 
                                 "second.best.HGT.taxon"], axis = 1)
rsoc_data.to_csv(r"C:\Users\Rhiannon\Documents\Year_4_Uni\project\python\rmag_GO.csv")

#%%
rsoc_missing = data.loc[list_of_bad]

for index, row in rsoc_missing.iterrows():
    print(index, row["go.numbers"], row["annotation"], row["best.HGT.taxon"], 
          row["second.best.HGT.taxon"], "\n")
    
rsoc_data = rsoc_missing.filter(["go.numbers", "annotation", "best.HGT.taxon", 
                                 "second.best.HGT.taxon"], axis = 1)
rsoc_data.to_csv(r"C:\Users\Rhiannon\Documents\Year_4_Uni\project\python\rsoc_GO.csv")
