# -*- coding: utf-8 -*-
"""
Created on Thu May 16 11:08:11 2019

@author: Rhiannon
"""

import pandas as pd
import csv
import os
from Bio import SeqIO
os.getcwd()
os.chdir("C:\\Users\\Rhiannon\\Documents\\Year_4_Uni\\project\\transcript_data")
#%%

transcriptIDlist = []
geneIDlist = []

rmag_genes_file = SeqIO.parse(open("HGT2.fasta"), "fasta")
for fasta in rmag_genes_file:
    name, sequence = fasta.description, str(fasta.seq)
    print(name)
    splitname = name.split(" ")
    transcriptIDlist.append(splitname[0])
    geneIDlist.append(splitname[4])
    
    #%%Aric
Aric_list = []
with open("Aric_discontig.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        rowstr = str(row["Name"])
        #print(rowstr)
        refSeq = row["% Of Ref Seq"]
        refSeq = refSeq[:-1]
        coverage = int(row["Min Sequence Length"])/int(row["Sequence Length"])
        print(coverage)
        if float(refSeq) > 40 and coverage > 0.5:
           Aric_list.append(rowstr)
        
Aric_pres = []
for ID in transcriptIDlist:
    if ID in Aric_list:
        Aric_pres.append("1")
    else:
        Aric_pres.append("0")

#%%Rmac
Rmac_list = []
print("Rmac")
with open("Rmac_discontig.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        rowstr = str(row["Name"])
        #print(rowstr)
        refSeq = row["% Of Ref Seq"]
        refSeq = refSeq[:-1]
        coverage = int(row["Min Sequence Length"])/int(row["Sequence Length"])
        print(coverage)
        if float(refSeq) > 40 and coverage > 0.5:
           Rmac_list.append(rowstr)
        
Rmac_pres = []
for ID in transcriptIDlist:
    if ID in Rmac_list:
        Rmac_pres.append("1")
    else:
        Rmac_pres.append("0")

        
#%%
trans_data = {
        "TranscriptID" : transcriptIDlist,
        "GeneID" : geneIDlist,
        "Aric" : Aric_pres,
        "Rmac" : Rmac_pres}
HGT_Aric_df = pd.DataFrame(trans_data, columns=['TranscriptID', 'GeneID', 
                               'Aric', 'Rmac'])
 
print("Transcript", HGT_Aric_df, sep='\n')

#%%
HGT_Aric_df['Aric'].value_counts()

#1    152
#0     86
HGT_Aric_df['Rmac'].value_counts()
#1    144
#0     94

#%%
dict_of_potential_transfer = {}
list_transID = []
for index, row in HGT_Aric_df.iterrows():
    #print(row["Aric"])
    if row["Aric"] == "0" and row["Rmac"]== "0":
        print( row["GeneID"], row["TranscriptID"])
        dict_of_potential_transfer[row["TranscriptID"]] = row["GeneID"]
        list_transID.append(row["TranscriptID"])
        
#%%
Aric_list_HGT = open("HGT_noAricRmac.fasta","w+")            
           
genes_file = SeqIO.parse(open("HGT2.fasta"), "fasta")
for fasta in genes_file:
    name, sequence = fasta.description, str(fasta.seq)
    print(name)
    splitname = name.split(" ")
    fasta_seq = str(">" + str(name) + "\n" + str(sequence) + "\n")  
    if splitname[0] in list_transID:
        Aric_list_HGT.write(fasta_seq)

        
Aric_list_HGT.close()