# Heatmap-78-4CC

This repository contains the main scripts to generate heatmaps and barplots of antimicrobial resistance genes (ARGs) from orthology analysis results. The analysis integrates bacterial isolate genomes, the CARD database, and ProteinOrtho outputs.

## Included Files

- `01-buf.py`: Python script for data processing.
- `02-cut.sh`: Bash script for file preparation.
- `heatmap-ortos2-SD.R`: R script to generate data visualizations.
- `R-env.yml`: R environment file listing package dependencies.
- `.gitignore`: Specifies untracked files and folders to exclude from the repository.

## Missing directories

Due to space restrictions, files are not included. 

  `01-Genomas_80_blactamasa_FAA/`: this directory contains the annotated genomes and CARD AA sequences.
  `02-proteinortho`: Generated with ProteinOrtho proteinortho6.pl -project=myproject 01-Genomas_80_blactamasa_FAA/*.faa > 02-proteinortho/myproject.proteinortho

## 🔧 Tools and Technologies

- R (e.g., `ggplot2`, `pheatmap`)
- Python
- Bash
- [ProteinOrtho](https://github.com/lschm/Proteinortho)
- [CARD Database](https://card.mcmaster.ca/)

## 🎯 Objective

To visualize patterns of antimicrobial resistance by comparing orthologous genes across bacterial genomes isolated from extreme environments.

## 🌐 Remote Repository

[https://github.com/mizas/heatmap-78-4CC](https://github.com/mizas/heatmap-78-4CC)

---

Feel free to expand this file with sections like "Usage", "Installation", or "Citing this work" depending on your needs.
