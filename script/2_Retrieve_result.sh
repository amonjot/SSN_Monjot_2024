#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 24/07/2024
#
mkdir -p result/Lagoon_output/
cp SSN_env/results/attributes/Metadata_Unicellular.attributes result/Lagoon_output/
for treshold in 80_65_1e-50 80_70_1e-50 80_75_1e-50 80_80_1e-50 80_85_1e-50 80_90_1e-50 80_95_1e-50 80_100_1e-50
do
    diamonond_out=$(echo -e "diamond_ssn_"$treshold".edges")
    cluster_out=$(echo -e "ssn_composantes_connexes_graph_"$treshold".csv")
    stat_out=$(echo -e "ssn_statistique_"$treshold".txt")
    mkdir result/Lagoon_output/$treshold
    cp SSN_env/results/$diamonond_out result/Lagoon_output/$treshold/
    cp SSN_env/results/$cluster_out result/Lagoon_output/$treshold/
    cp SSN_env/results/$stat_out result/Lagoon_output/$treshold/
done
