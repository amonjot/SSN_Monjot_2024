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
# Join TranscriptID and ProteinID table
# Import modules
import re
import os
import math
from collections import defaultdict

#Taxonomyold2Taxonomynew
# Declare Dictionnary
transID2protID=defaultdict(str)
# Read
with open ("rawdata/proteins_from_Unigenes_CEQ.fa","r") as f1:
    for li in f1 :
        li=li.rstrip("\n")
        if li.startswith(">"):
            ProteinID=li.lstrip(">")
            TranscriptsID=ProteinID.split("-")[0]
            transID2protID[TranscriptsID]=ProteinID

# TranscriptID2Metadata
# Declare Dictionnary
transID2METADATA=defaultdict(str)
# Read
with open ("rawdata/out_RESUME.txt","r") as f2:
    for li in f2 :
        li=li.rstrip("\n")
        ID_TRANSCRIPT=li.split("\t")[0]
        METADATA="\t".join(li.split("\t")[1:])
        transID2METADATA[ID_TRANSCRIPT]=METADATA
#
# Transcript2ko
# Declare Dictionnary
protID2ko_dict=defaultdict(str) # Dictionnary associating transcripts to taxonomy
# Read
with open ("rawdata/annotation_ko.tsv","r") as f7:
    for li in f7 :
        li=li.rstrip("\n")
        if li.startswith("CEQ"):
            ID_prot="-".join(li.split("\t")[0:2])
            ko=li.split("\t")[2]
            protID2ko_dict[ID_prot]=ko
#
#
# Create result table
for TranscriptsID in transID2protID:
    print(transID2protID[TranscriptsID]+"\t"+TranscriptsID+"\t"+transID2METADATA[TranscriptsID]+"\t"+protID2ko_dict[transID2protID[TranscriptsID]])
#
