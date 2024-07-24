#!/usr/bin/env bash

fasta=${1}
db=${2}
alignment=${3}
th=${4}

diamond makedb --in $fasta --threads $th --db $db

diamond blastp -d $db \
	       -q $fasta \
	       -o $alignment \
	       -e 1e-5 \
	       --sensitive \
           --threads $th \
	       --outfmt 6 qseqid sseqid pident ppos length mismatch gapopen qstart qend sstart send evalue bitscore
