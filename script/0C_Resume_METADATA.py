#!/usr/bin/env python
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 14/04/2023
#
# Resume metadata table
# Import modules
import re
import os
import math
from collections import defaultdict

#Taxonomyold2Taxonomynew
# Declare Dictionnary
taxold2new_dict=defaultdict(str) # Dictionnary associating transcripts to taxonomy
# Read
with open ("rawdata/Taxonomy_correspondance.txt","r") as f2:
    for li in f2 :
        li=li.rstrip("\n")
        OLD=li.split("\t")[0]
        NEW=li.split("\t")[1]
        taxold2new_dict[OLD]=NEW
#
#Taxonomynew2TRT
#Declare Dictionnary
taxnew2TRT_dict=defaultdict(str)
# Read
with open ("rawdata/Table_S1.tsv","r") as f3:
    for li in f3 :
        li=li.rstrip("\n")
        TAXO=li.split("\t")[0]
        SPCS=li.split("\t")[2]
        TRT1=li.split("\t")[3]
        TRT2=li.split("\t")[4]
        TYP1=li.split("\t")[5]
        TYP2=li.split("\t")[6]
        taxnew2TRT_dict[TAXO]=SPCS+"\t"+TRT1+"\t"+TRT2+"\t"+TYP1+"\t"+TYP2
#
# Transcript2Taxonomy
# Declare Dictionnary
trans2Tax_dict=defaultdict(str) # Dictionnary associating transcripts to taxonomy
# Read
with open ("rawdata/table_taxonomy.perUnigene.allUnigenes.tsv","r") as f1:
    for li in f1 :
        li=li.rstrip("\n")
        if li.startswith("CEQ"):
            TAXO=li.split("\t")[1]
            ID_TRANSCRIPT=li.split("\t")[0]
            trans2Tax_dict[ID_TRANSCRIPT]=TAXO
#
#
#
# Create result table
for ID_TRANSCRIPT in trans2Tax_dict:
    print(ID_TRANSCRIPT+"\t"+trans2Tax_dict[ID_TRANSCRIPT]+"\t"+taxold2new_dict[trans2Tax_dict[ID_TRANSCRIPT]]+"\t"+taxnew2TRT_dict[taxold2new_dict[trans2Tax_dict[ID_TRANSCRIPT]]])
#
#
#
