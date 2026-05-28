# Documentation of the article: Descovering the genomic landscape of liver hepatocarcinoma across mutational signatures and driver genes

En este documento se presenta la informacion, referencias, figuras intermedias y codigos utilizados para escribir el articulo correspondiente

- **Authors:**
  - Román Cervantes Levario

  -  Jair Emilianto Contreras Rivera
  
  - Leonardo Daniel González López

- **Date:** 20/05/2026

- **Semester:** 4th semester

## Abstract

Liver hepatocellular carcinoma (LIHC) is the most common primary liver cancer, accounting for approximately 90% of all cases. The accumulation of somatic driver mutations and characteristic mutational signatures are key mechanisms underlying its development. However, the genomic landscape of LIHC remains incompletely characterized.

Here, we performed a genomic analysis of 358 LIHC tumor samples using a MAF file containing 60,145 mutations obtained by whole-exome sequencing. Driver gene and hotspot detection was performed using dndscv and mutational signature analysis using the SigProfiler framework mapped to COSMIC reference signatures.

We identified 15 driver gene candidates, including the oncogenes *CTNNB1* and *NFE2L2*, and 9 tumor suppressors, with *TP53* being the most frequently mutated. Hotspot analysis consistently identified 7 missense hotspots in *CTNNB1* exon 3 and 2 hotspots in *TP53* affecting the DNA-binding domain. Mutational signature analysis revealed signatures 5, 22a and 40a as the most active, with signature 22a, associated with aristolochic acid exposure, showing a distinct distribution among driver gene mutated samples.

These findings highlight the value of integrating algorithms for driver genes detection and mutational signatures analysis with data platforms to characterize tumor biology, identify potential therapeutic targets, and understand the molecular mechanisms of hepatocellular carcinoma.

## R packages
```
# Driver genes
library(devtools)
library(dndscv)

# Mutational signatures
library(tidyverse)
library(ggplot2)
library(reticulate)
library(devtools)
library(SigProfilerMatrixGeneratorR)
library(SigProfilerPlottingR)
library(SigProfilerExtractorR) 
library(SigProfilerAssignmentR)
```
