#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 25/09/2023
#
BEFORE=$SECONDS
# R CStack
ulimit -s unlimited
# Activate conda environment
PATHCONDA=$(conda info | grep -i 'base environment' | awk -F" " '{print $4}')
source $PATHCONDA'/etc/profile.d/conda.sh'
conda activate REnv_Monjot_2024A
if [ $(echo $1 | grep "." |wc -l ) == 1 ]
then
OUTPUT=$1
else
echo 'Enter result file path : '
read OUTPUT
fi

echo "Output: "$OUTPUT


# set script directory
mkdir Monjot_etal_2024
mkdir Monjot_etal_2024/Figures
mkdir Monjot_etal_2024/Figures_Supp
#Fig 1
##A
cp result/Lagoon_output/$OUTPUT/Analysis/Stat_In_CC.svg Monjot_etal_2024/Figures/Figure_1A.svg
##B
cp result/Lagoon_output/$OUTPUT/Analysis/Stat_bis_In_CC.svg Monjot_etal_2024/Figures/Figure_1B.svg
#Fig 2
cp result/Lagoon_output/$OUTPUT/Analysis/Grp_Intersec_in_CC_UPSET_data.svg Monjot_etal_2024/Figures/Figure_2.svg
#Fig 3
cp result/Lagoon_output/$OUTPUT/Analysis/CCA_TRT.svg Monjot_etal_2024/Figures/Figure_3.svg
#Fig 4
cp result/Lagoon_output/$OUTPUT/Analysis/taxo_compo_75.svg Monjot_etal_2024/Figures/Figure_4.svg
#Fig 6
##A
cp result/Lagoon_output/$OUTPUT/Igraph/CC_containing75PARA_CC_PARA_obligate.svg Monjot_etal_2024/Figures/Figure_6A.svg
##B
cp result/Lagoon_output/$OUTPUT/Igraph/CC_containing75PARA_CC_PARA_obligate_KO_50most_abundant.svg Monjot_etal_2024/Figures/Figure_6B.svg
#Fig 7
cp result/Lagoon_output/$OUTPUT/Analysis/CCA_TRT_Barycentre_obligatePARA.svg Monjot_etal_2024/Figures/Figure_7.svg
#
#Fig S1
cp result/Lagoon_output/Parameter_analyses.svg Monjot_etal_2024/Figures_Supp/Figure_S1.svg
#Fig S2
cp result/Lagoon_output/$OUTPUT/Analysis/Number_Prot_by_CC.svg Monjot_etal_2024/Figures_Supp/Figure_S2.svg
#Fig S3
cp result/Lagoon_output/$OUTPUT/Analysis/CCA_TRT_eig.svg Monjot_etal_2024/Figures_Supp/Figure_S3.svg
#Fig S4
cp result/Lagoon_output/$OUTPUT/Analysis/KO_compo_75.svg Monjot_etal_2024/Figures_Supp/Figure_S4.svg
#Fig S5
cp result/Lagoon_output/$OUTPUT/Analysis/KO_compo_special_category.svg Monjot_etal_2024/Figures_Supp/Figure_S5.svg
#Fig S6
cp result/Lagoon_output/$OUTPUT/Analysis/Stat_CC_75percent_obligate.svg Monjot_etal_2024/Figures_Supp/Figure_S6.svg
#Fig S7
cp result/Lagoon_output/$OUTPUT/Igraph/CC_containing75PARA_CC_PARA_obligate_KO_correlated.svg Monjot_etal_2024/Figures_Supp/Figure_S7.svg


ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "R process finished and takes "$ELAPSED" minutes !"
