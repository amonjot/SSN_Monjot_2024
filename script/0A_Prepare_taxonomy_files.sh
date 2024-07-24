#!/bin/bash
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
# Ordered taxon
## line to line reading
OLDIFS=$IFS
IFS=$'\n'
##
cat rawdata/table_taxonomy.perUnigene.allUnigenes.tsv | grep -v "Unigene\tValues" | awk -F'\t' '{print $2}' | sort | uniq > rawdata/taxon.txt
## Start
for affiliation in $(cat rawdata/taxon.txt)
do
### set variable
Domain=$(echo $affiliation | awk -F';d_' '{ print $2 }' | awk -F';' '{ print $1 }')
Division=$(echo $affiliation | grep ";d_Eukaryota;-_" | awk -F';d_Eukaryota;-_' '{ print $2 }' | awk -F';' '{ print $1 }')
Kingdom=$(echo $affiliation | grep ";k_" | awk -F';k_' '{ print $2 }' | awk -F';' '{ print $1 }')
Phylum=$(echo $affiliation | awk -F';p_' '{ print $2 }' | awk -F';' '{ print $1 }')
Class=$(echo $affiliation | awk -F';c_' '{ print $2 }' | awk -F';' '{ print $1 }')
Order=$(echo $affiliation | awk -F';o_' '{ print $2 }' | awk -F';' '{ print $1 }')
Family=$(echo $affiliation | awk -F';f_' '{ print $2 }' | awk -F';' '{ print $1 }')
Genus=$(echo $affiliation | awk -F';g_' '{ print $2 }' | awk -F';' '{ print $1 }')
Species=$(echo $affiliation | awk -F';s_' '{ print $2 }' | awk -F';' '{ print $1 }')
### echo results
if [ $(echo $Division | wc -c) -gt 2 ]
then
	if [ $(echo $Kingdom | wc -c) -gt 2 ]
	then
	echo -e $affiliation'\t'$Domain'/d_'$Division'/k_'$Kingdom'/p_'$Phylum'/c_'$Class'/o_'$Order'/f_'$Family'/g_'$Genus >> rawdata/Taxonomy_correspondance.txt
	elif [ $(echo $Kingdom | wc -c) -eq 1 ]
	then
	echo -e $affiliation'\t'$Domain'/d_'$Division'/p_'$Phylum'/c_'$Class'/o_'$Order'/f_'$Family'/g_'$Genus >> rawdata/Taxonomy_correspondance.txt
	fi
elif [ $(echo $Division | wc -c) -eq 1 ]
then
echo -e $affiliation'\t'$Domain'/p_'$Phylum'/c_'$Class'/o_'$Order'/f_'$Family'/g_'$Genus >> rawdata/Taxonomy_correspondance.txt
fi
done
#
IFS=$OLDIFS

# Sort only uniq taxonomy
cat rawdata/Taxonomy_correspondance.txt | awk -F'\t' '{ print $2}' | sort | uniq > rawdata/unique.csv

# detect last taxonomy
## line to line reading
OLDIFS=$IFS
IFS=$'\n'
## Start
# Define $Last
for affiliation in $(cat rawdata/unique.csv)
do
#echo $affiliation
i=$(($(echo $affiliation | grep -o "/" | wc -l)+1))
Last_temp=$(echo $affiliation | awk -v temp=$i -F'/' '{print $temp}')
while [ $(echo $Last_temp | grep "_$" | wc -l) -eq 1 ]
do
i=$(($i-1))
Last_temp=$(echo $affiliation | awk -v temp=$i -F'/' '{print $temp}')
done
Last=$(echo $Last_temp | awk -F'_' '{print $2}')
# print result
echo -e $affiliation"\t"$Last >> rawdata/Tax_Table.tsv
done
#
IFS=$OLDIFS

## Clean directory
rm rawdata/unique.csv

# Search Species
## line to line reading
OLDIFS=$IFS
IFS=$'\n'
#
for i in $(cat rawdata/Tax_Table.tsv)
do
Last=$(echo $i | awk -F"\t" '{ print $2}')
Genus=$(echo $i | awk -F"\t" '{ print $1}' | awk -F"/g_" '{ print $2}')
#echo "Last:"$Last
#echo "Genus:"$Genus
speciesf=""
#
if [ $(echo $Genus | wc -c) -gt 1 ]
then
for species in $(cat rawdata/taxon.txt | grep "g_$Genus" | awk -F";s_" '{ print $2 }' | cut -d";" -f1 | sort | uniq)
do
speciesf=$(echo $speciesf" + "$species)
done
speciesf=$(echo $speciesf | cut -c 3-)
fi
#
echo -e $i"\t"$(echo $speciesf | sed 's/^ //g') >> rawdata/Table_assoc_Temp.tsv
done

IFS=$OLDIFS

## Clean directory
rm rawdata/taxon.txt
rm rawdata/Tax_Table.tsv
