
# Overview

This is the accompanying code repository for the paper Schropp et al. (under revision) - "The Impact of Regular Sauerkraut Consumption on the Human Microbiota: a Cross-Over Intervention Trial".

The code contains the data and R code for all relevant analyses in the paper, including the supplementary material.
Data and code are published under the CC BY 4.0 license.


# Folder structure

- `code (helpers)`: Helper functions for the main code, in form of an R package. The package can be locally installed by calling `devtools::install("code (helpers)/sauerkrautTaxonomyBuddy")`.
- `code (main)`: Main R code files (using Quarto reports) reproducing all analyses and figures shown in the paper.
- `data`: Data folder (see descriptions below)
- `figures`: Figure files produced by the code.


# Data description

- `1a - group info.Rdata`: Information to which intervention group (A or B, receiving the fresh or pasteurized sauerkraut intervention first, respectively) each participant belongs.
- `1b - participant data.Rdata`: Several datasets containing general participant information like which intervention was received at which timepoint and general information like BMI, weight, etc.
- `1c - sample IDs.xlsx`: General base information on observations, including some information on duplicate measurements (which were measured in duplicate for estimating the measuring reliability)
- `1d - food frequency questionnaire.csv`: Food questionnaire responses
- `1e - food frequency questionnaire (nutrients).csv`: Nutrient data based on the food questionnaire responses
- `1f - opinion questionnaire.xlsx`: Opinion questionnaire responses, with explanations of the items in file `1f - lookup for opinion questionnaire.xlsx`.
- `1g - blood markers.xlsx`: Blood marker measurements (glucose, CRP, etc.)
- `1h - body measurements.xlsx`: Body measurements (BMI, blood pressure, etc.)
- `1i - symptom diaries.xlsx`: Symptom diary data
- `2a - stool taxonomy (absolute abundances).tsv`: Taxonomic stool measurements (absolute abundances)
- `2b - stool taxonomy (relative abundances).tsv`: Taxonomic stool measurements (relative abundances)
- `2c - meta information blood and stool samples.xlsx`: Meta information on blood and stool samples and duplicate measurements (measured for estimating measuring reliability)
- `2d - meta information sauerkraut samples.xlsx`: Some (quite unstructured) general information on the taxonomic samples which were no human stool samples, but taxonomic samples taken from some sauerkraut jars
- `3a - KEGG genes.txt`: KEGG KO data on genes in the stool samples
- `3b - KEGG pathways.txt`: KEGG data on pathways in the stool samples
- `4a - blood data part 1.Rdata`: Further blood measurements (TNFR2, RAGE, etc.)
- `4b - blood data part 2.Rdata`: Further blood measurements (glucose, sCRP, etc.)
- `4c - blood data part 1 including duplicates.xlsx`: Blood measurements (TNFR2, RAGE, etc.), also including the duplicate samples (used for estimating measuring reliability)
- `4d - blood data part 2 including duplicates.xlsx`: Blood measurements (glucose, sCRP, etc.), also including the duplicate samples (used for estimating measuring reliability)
- `5a - blood SCFAs.xlsx`: Blood SCFA measurements
- `5b - stool SCFAs.xlsx`: Stool SCFA measurements
