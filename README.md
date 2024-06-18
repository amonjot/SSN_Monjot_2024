# **Metatranscriptomes-based sequence similarity networks uncover genetic signatures within parasitic freshwater microbial eukaryotes**

**This is a workflow to reproduce analysis conduced in Monjot et al., 2024**

First, clone github repository: `git clone https://github.com/amonjot/SSN_Monjot_2024.git`

Second, define current directory: `cd SSN_Monjot_2024`

                                        *******************************
                                            Directory organization
                                        *******************************
    Microstore_Analysis_Monjot_et_al._2023 
    |-> rawdata (sub-directory for trophic modes table and metadata)
        |-> ko_to_hierarchy.txt (KO id definition table from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg)
        |-> Table_S1.tsv (own in-house trophic modes database (Monjot et al., 2024))
        |-> annotation_ko.tsv (The output of the KO annotation step using koFamScan) [https://doi.org/10.5281/zenodo.11972386]
        |-> proteins_from_Unigenes_CEQ.fa (The proteins fasta file provided by the proteins prediction step after filtration) [https://doi.org/10.5281/zenodo.11972386]
        |-> table_taxonomy.perUnigene.allUnigenes.tsv (The output of the taxonomic affiliation against MetaEuk) [https://doi.org/10.5281/zenodo.11972386]
    |-> script (sub-directory for minor scripts)
        |-> 1_Choose_Treshold.R (script to test all identity treshold and determine functional homogeneity)
        |-> 2_Igraph (script to process CCs graph)
        |-> 0A_parse_ko_hits.py (script to parse ko annotation using the first hit) 
        |-> 0B_Prepare_taxonomy_files.sh (script to Prepare the trophic mode database model)
        |-> 0C_Resume_METADATA.py (script to resume data related to transcripts)
        |-> 0D_Cross-TranscriptID-ProteinID.py (script to link transcripts data and proteins ID)
        |-> environment_REnv_Monjot_2024A.yml (conda environment)
        |-> REnv_Monjot_2023A_packages (local repository to R packages installation)
        |-> Preprocess_setup.sh (script to launch conda environment setup)
        |-> 4_Retrieve_Figures.sh (script to retrieve published figures from result directory)

    
## 1. Conda environment and dependencies installation

### Conda and dependencies
Install miniconda following the standard procedure (https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html)

Then, install conda environment with the following script: 

    bash Preprocess_setup.sh

* script: Preprocess_setup.sh
* This takes 1 min on [Dual CPU] Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

This installs the following tools:

    * KronaTools v2.8
    * cutadapt v4.1
    * r-base v4.2.2
    * imagemagick v7.1.0_52
    * python v3.9.15
    * perl v5.32.1
    * cmake v3.22.1
    * parallel v20230722
    * numerous libraries (details in script/environment_REnv_Monjot_2023A.yml)

### R installation (+ dependencies) if necessary

Install R packages: 

    * ggplot2           * ggstatsplot
    * plyr              * stringr
    * dplyr             * treemapify
    * tidyr             * elementalist
    * cowplot           * gplots
    * paletteer         * igraph
    * reshape2          * ggraph
    * varhandle         * tidygraph
    * ggrepel           * grid
    * ggpubr            * FactoMineR
    * ggsci             * factoextra
    * scales            * rstatix
    * hrbrthemes        * ggh4x
    * svglite           * fmsb
    * ggupset


## 2. Metatranscriptomics data - filtering and annotation

### Raw data: availability, assembly and cleanning

The raw data are available in the public databases under the umbrella of the _BioProject PRJEB61515_.    
There are 32 samples, each of them assembled one-by-one
with _oases_ and were clustered all together with _CD-HIT-EST_ with the thresholds
identity &gt;95% over &gt;90% of the length of the smallest sequence. Moreover,
transcripts longer than 50kb were discarded, this resulted in 10.359.104
representative assembled transcripts, called herafter _Unigenes_. The protocol
is described in grater details in the work from
[Carradec _et al._ 2018](https://doi.org/10.1038/s41467-017-02342-1).

We then removed human contamination from the _Unigenes_. _Unigenes_ were aligned
to the Human genome assembly _GRCh38.p13_ with _minimap2 v2.14-r894-dirty_ and
the default parameters. All _Unigenes_ with a hit against the genome weere
considered as contaminant as thus discarded for future analysis. Here is an example
of code to obtain the list of contaminants:

```bash
minimap2 -t 12 GRCh38.p13.genome.fa.gz Unigenes.fa | cut -f 1 | sort | \
    uniq >list_unigene_human_contaminant
```

A total of 32.218 _Unigenes_ were discarded.

### Protein prediction

We used [_TransDecoder v5.5.0_](https://github.com/TransDecoder/TransDecoder/wiki)
to predict coding sequences present on the _Unigenes_. The minimal protein length
was set to 70 amino-acids as we observed _Unigenes_ without protein with the
default parameters:

```bash
# Extract long Open-Reading Frames 
TransDecoder.LongOrfs -m 70 --output_dir out_transDecoder_long -t unigenes.fa

# Predict the likely coding regions
TransDecoder.Predict --output_dir out_transDecoder_pred -t unigenes.fa
```

Proteins have been check against the [_AntiFam_](https://doi.org/10.1093/database/bas003)
database and _HMMER v3.3.2_ using the profiles' score cutoff. Spurious proteins
were then discarded:

```bash
# Resources:
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz
tar -zxf Antifam.tar.gz

# Run the comparison: only positive hits
hmmsearch --cut_ga --noali --tblout antifam_search.tsv AntiFam.hmm proteins.fa
```

### _Unigenes_ taxonomic affiliation

We used _MetaEuk_ version _commit 57b63975a942fbea328d8ea39f620d6886958eca_.
The reference for the taxonomic affiliation is the database provided by _MetaEuk_
authors, available
[here](https://wwwuser.gwdguser.de/~compbiol/metaeuk/2020_TAX_DB/). This web-page
proposes a link to download the data as well as a description of the origin of data
composing this dataset. Beware the database is a 20 GB _tar.gz_ archive that takes
up to **200 GB** of disk-space once uncompressed.

```bash
MetaEukTaxoDB=MMETSP_zenodo_3247846_uniclust90_2018_08_seed_valid_taxids

# Create the 'MMSeqs' database
metaeuk createdb unigenes.cleaned.fa UnigeneDB

# Search the taxonomy for each protein
metaeuk taxonomy UnigeneDB $MetaEukTaxoDB Unigene_taxoDB tmp \
    --majority 0.5 --tax-lineage 1 --lca-mode 2 --max-seqs 100 -e 0.00001 \
    -s 6 --max-accept 100

# Get a Krona-like plot
metaeuk taxonomyreport $MetaEukTaxoDB Unigene_taxoDB \
    Unigene_report.html --report-mode 1

# Get a tsv
metaeuk createtsv UnigeneDB Unigene_taxoDB \
    Unigene_taxonomy_result.tsv
```

Then we associated the taxonomy of the protein to its corresponding _Unigene_.
In cases where several proteins were available for a single _Unigene_, 


#### Clean contaminant from 

From the taxonomic information, _Unigenes_ with proteins affiliated to _Bacteria_,
_Archaea_ or _Viruses_ were removed ==> 249 878 Unigenes discarded

Metazoans too, same way with the taxonomic affiliation of proteins associated to
Unigenes: 150.368 _Unigenes_ removed

single protein ==> report the affiliation to the _Unigene_  
several proteins ==> report the LCA

### Proteins annotations

#### KEGG _KO_

Proteins have been annotated with the KEGG's _KO_ through the tool
[_koFamScan v1.3.0_](https://www.genome.jp/tools/kofamkoala/) using the _KO_
HMM profiles release "_2022-01-03_", available
[here](https://www.genome.jp/ftp/db/kofam/).

```bash
# Get data
wget https://www.genome.jp/ftp/db/kofam/archives/2022-01-03/ko_list.gz
gunzip ko_list.gz
wget https://www.genome.jp/ftp/db/kofam/archives/2022-02-01/profiles.tar.gz
tar -zxf profiles.tar.gz

# Run
exec_annotation -o results.koFamScan.tsv --format detail-tsv --ko-list ko_list\
    --profile profiles proteins.fa
```

Then, we sent the results to the _Python3_ script
`script/parse_ko_hits.py`:

```bash
python3 script/parse_ko_hits.py --input results.koFamScan.tsv \
    --output result/results.koFamScan.parsed.tsv
```

This script parses the results in this order of preference:

1. Keep hit if tagged as significant by _KoFamScan_, the ones with a `*` in the
first field. This type of result is tagged with "_significant_" in the result
2. If the current protein has **no** significant hit, keep the best hit if its
_e-value_ is &le; 1.e-5.

Results of annotation process are available [here](https://doi.org/10.5281/zenodo.11972386).
Download these files and place it in the rawdata directory for use.

## 3. Trophic mode database and SSN inputs

### Trophic mode database

We generate a table ready to be filled corresponding to the model of the trophic mode database : 

```bash
bash script/0B_Prepare_taxonomy_files.sh
```

Then, we have to complete the resulted table Table_assoc_Temp.tsv using Bibliography to create the final trophic mode database.
This table correspond to the Supplementary file *Table_S1.tsv* and is provided in the rawdata directory.

This table was developed as part of a metatranscriptomic study (Monjot et al., 2023) using bibliography (362 scientific articles), 
to link, when possible, proteins taxonomic affiliation (i.e. 1 052 distinct affiliations), to trophic modes.

Taxonomic ranks are specified with the following code: d_(Division), k_(Kingdom), p_(Phylum), c_(Class), o_(Order), f_(Family) and g_(Genus). 
Trophic modes were identified using the following codes: PHOTO: photo-osmo-mixotrophs, MIXO: photo-osmo-phago-mixotrophs, HET: heterotrophs, 
SAP: saprotrophs and PARA: parasites (facultative and obligate).

### Metadata files
 
Cross all proteins informations: 

 ```bash
python3 script/0C_Resume_METADATA.py > rawdata/out_RESUME.txt
python3 script/0D_Cross-TranscriptID-ProteinID.py > rawdata/out-RESUME-PROTID.txt
```

Only informations characterizing proteins affiliated to microbial eukaryotes should be retained (by removing those affiliated to pluri- or multi-cellular organisms):

```bash
echo -e "ID_Protein"\\t"ID_transcript"\\t"ncbi_taxonomy"\\t"Taxonomy"\\t"SpeciesInDataset"\\t"TRT_1"\\t"TRT_2"\\t"TYP_1"\\t"TYP_2"\\t"ko" > Metadata_Unicellular.txt
cat rawdata/out-RESUME-PROTID.txt | grep -v "Pluri- or multi-cellular" >> Metadata_Unicellular.txt
```

### Fasta file

From the inital fasta file, we remove all sequences affiliated to pluri- or multi-cellular organisms:

```bash
cat rawdata/out-RESUME-PROTID.txt | grep -v "Pluri- or multi-cellular" | cut -f1 > rawdata/trsc_Unicellular.txt
seqtk subseq rawdata/proteins_from_Unigenes_CEQ.fa rawdata/trsc_Unicellular.txt > proteins_from_Unigenes_CEQ_Unicellular.fasta
```

## 4. SSN analysis



## 5. Retrieve article figures

To retrieve article figures, run following script:

    bash script/4_Retrieve_Figures.sh 80_80_1e-50

* script: 4_Retrieve_Figures.sh
* argument 1: 80_80_1e-50
* This takes just a few seconds on [Dual CPU] Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

The resulting figures can be found in *Monjot_etal_2024/* directory.
