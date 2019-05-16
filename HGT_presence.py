# -*- coding: utf-8 -*-
"""
Created on Tue May 14 16:40:06 2019

@author: Rhiannon
"""
#%%
import pandas as pd
import csv
import os
from Bio import SeqIO
os.getcwd()
os.chdir("C:\\Users\\Rhiannon\\Documents\\Year_4_Uni\\project\\transcript_data")
#%%

transcriptIDlist = []
geneIDlist = []

rmag_genes_file = SeqIO.parse(open("rmag_extract.txt"), "fasta")
for fasta in rmag_genes_file:
    name, sequence = fasta.description, str(fasta.seq)
    print(name)
    splitname = name.split(" ")
    transcriptIDlist.append(splitname[0])
    geneIDlist.append(splitname[4])
    
rmag_genes_file = SeqIO.parse(open("rsoc_extract.txt"), "fasta")
for fasta in rmag_genes_file:
    name, sequence = fasta.description, str(fasta.seq)
    print(name)
    splitname = name.split(" ")
    transcriptIDlist.append(splitname[0])
    geneIDlist.append(splitname[4])
        

#%%Aric
Aric_list = []
with open("Aric_hits.csv") as myFile:
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

#%%Avag
        
Avag_list = []
with open("Avag_hits.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        rowstr = str(row["Name"])
        print(rowstr)
        refSeq = row["% Of Ref Seq"]
        refSeq = refSeq[:-1]
        coverage = int(row["Min Sequence Length"])/int(row["Sequence Length"])
        print(coverage)
        if float(refSeq) > 40 and coverage > 0.5:  
           Avag_list.append(rowstr)
        
Avag_pres = []
for ID in transcriptIDlist:
    if ID in Avag_list:
        Avag_pres.append("1")
    else:
        Avag_pres.append("0")

#%%Rmac

Rmac_list = []
with open("Rmac_hits.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        rowstr = str(row["Name"])
        print(rowstr)
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

#%%Rmag
Rmag_list = []
with open("Rmag_hits.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        rowstr = str(row["Name"])
        print(rowstr)
        refSeq = row["% Of Ref Seq"]
        refSeq = refSeq[:-1]
        coverage = int(row["Min Sequence Length"])/int(row["Sequence Length"])
        print(coverage)
        if float(refSeq) > 40 and coverage > 0.5:
           Rmag_list.append(rowstr)
        
Rmag_pres = []
for ID in transcriptIDlist:
    if ID in Rmag_list:
        Rmag_pres.append("1")
    else:
        Rmag_pres.append("0")
        
#%%Rsoc
Rsoc_list = []
with open("Rsoc_hits.csv") as myFile:
    reader = csv.DictReader(myFile)
    for row in reader:
        rowstr = str(row["Name"])
        print(rowstr)
        refSeq = row["% Of Ref Seq"]
        refSeq = refSeq[:-1]
        coverage = int(row["Min Sequence Length"])/int(row["Sequence Length"])
        print(coverage)
        if float(refSeq) > 40 and coverage > 0.5:
           Rsoc_list.append(rowstr)
        
Rsoc_pres = []
for ID in transcriptIDlist:
    if ID in Rsoc_list:
        Rsoc_pres.append("1")
    else:
        Rsoc_pres.append("0")

#%%
        
trans_data = {
        "TranscriptID" : transcriptIDlist,
        "GeneID" : geneIDlist,
        "Aric" : Aric_pres,
        "Avag" : Avag_pres,
        "Rmac" : Rmac_pres,
        "Rmag" : Rmag_pres,
        "Rsoc" : Rsoc_pres}
HGT_df = pd.DataFrame(trans_data, columns=['TranscriptID', 'GeneID', 
                               'Aric', 'Avag', 'Rmac', 'Rmag', 'Rsoc'])
 
print("Transcript", HGT_df, sep='\n')

#example
#dfObj = dfObj.append({'User_ID': 23, 'UserName': 'Riti', 'Action': 'Login'}, ignore_index=True)
#%%

HGT_df.to_csv(r"C:\Users\Rhiannon\Documents\Year_4_Uni\project\transcript_data\HGT_dataframe3.csv")

#%%
dict_of_potential_transfer = {}
list_transID = []
for index, row in HGT_df.iterrows():
    #print(row["Aric"])
    if row["Aric"] == "0" and row["Avag"]== "1" and row["Rmac"]== "0":
        print( row["GeneID"], row["TranscriptID"])
        dict_of_potential_transfer[row["TranscriptID"]] = row["GeneID"]
        list_transID.append(row["TranscriptID"])
        
#%%
full_list_HGT = open("HGT2.fasta","w+")            
           
rmag_genes_file = SeqIO.parse(open("rmag_extract.txt"), "fasta")
for fasta in rmag_genes_file:
    name, sequence = fasta.description, str(fasta.seq)
    print(name)
    splitname = name.split(" ")
    fasta_seq = str(">" + str(name) + "\n" + str(sequence) + "\n")  
    if splitname[0] in list_transID:
        full_list_HGT.write(fasta_seq)
        
    
rmag_genes_file = SeqIO.parse(open("rsoc_extract.txt"), "fasta")
for fasta in rmag_genes_file:
    name, sequence = fasta.description, str(fasta.seq)
    print(name)
    splitname = name.split(" ")
    fasta_seq = str(">" + str(name) + "\n" + str(sequence) + "\n")    
    if splitname[0] in list_transID:
        full_list_HGT.write(fasta_seq)
        
full_list_HGT.close()
                   




