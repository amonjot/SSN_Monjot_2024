# scripts4arthur

## Utilisation

Ce pipeline a été testé avec Ubuntu 22.04.

## Environnement

Afin de pouvoir exécuter chacune des étapes, il est nécessaire de construire un environnement conda avec environment.yaml. Ce fichier conteient tous les outils, langages et modules (python) utilisé  dans ce pipeline.

```bash
conda env create -f environment.yml

conda activate environment
```

## Utilisation

Ce pipeline utilise Snakemake. 
Snakemake est installer dans l'environnment Conda.

```bash
snakemake -c 10
```

