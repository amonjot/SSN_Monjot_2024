#!/bin/bash
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 18/06/2024
#
#
#
BEFORE=$SECONDS
PATHCONDA=$(conda info | grep -i 'base environment' | awk -F" " '{print $4}')
source $PATHCONDA'/etc/profile.d/conda.sh'
conda env create -f environment_REnv_Monjot_2024A.yml ; conda activate REnv_Monjot_2024A
#
ELAPSED=$((($SECONDS-$BEFORE)/60))
echo "Preprocess setup finished and takes "$ELAPSED" minutes !"
