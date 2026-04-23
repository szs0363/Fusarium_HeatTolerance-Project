# Fusarium_HeatTolerance-Project
FSSC Clade 3 genome analysis to identify unique genes to heat tolerant species.

## Overview

This project performs a comparative genomics analysis of **12 genomes from the *Fusarium solani species complex (FSSC) clade 3*** to identify orthogroups associated with **heat tolerance**.

Genomes were processed, annotated, and analyzed using a reproducible pipeline involving **Funannotate**, **OrthoFinder**, and downstream analysis in **R**.

---

## Data Source

* Genomes were obtained from **MycoCosm (JGI Genome Portal)**
* Dataset includes:

  * 12 FSSC Clade 3 genomes
  * Classified into:

    * Heat-tolerant species
    * Non-heat-tolerant species (based on metadata)

---

## Workflow Summary

### 1. Genome Processing

* Raw genomes were:

  * Cleaned
  * Sorted
  * Masked
* Tool used: **Funannotate (clean, sort, mask)**

---

### 2. Gene Prediction

* Protein-coding genes were predicted using:

  * `funannotate predict`
* Executed on:

  * High Performance Computing (HPC) cluster

---

### 3. Orthology Inference

* Input: Predicted protein FASTA files
* Tool: **OrthoFinder**
* Output:

  * Orthogroups
  * Gene count matrices
  * Phylogenetic inference files

---

### 4. Comparative Analysis

* Orthogroups were analyzed using metadata to identify:

  * Orthogroups unique to **heat-tolerant species**
  * Shared vs unique gene families

---

### 5. Visualization & Statistical Analysis (R)

R packages used:

* `tidyverse`
* `ggplot2`
* `pheatmap`
* `RColorBrewer`
* `UpSetR`

Analyses include:

* Presence/absence matrices
* Heatmaps of orthogroup distribution
* Venn/UpSet plots for shared gene families

---


## Outputs

Key outputs include:

* Orthogroup assignments
* Heat tolerance–specific gene families
* Heatmaps and clustering
* UpSet visualizations

---

## Key Objective

To identify **candidate orthogroups associated with heat tolerance** in FSSC clade 3 species using comparative genomics.

---

## Reproducibility Notes

* All scripts used in the pipeline are included in `/scripts`
* R scripts are modular and reproducible
* Metadata file is essential for grouping species

---

## Future Improvements

* Functional annotation of candidate orthogroups (e.g., InterProScan)
* Integration with differential expression data
* Identification of accessory genome contributions

---

[Link to my github](https://github.com/szs0363/Fusarium_HeatTolerance-Project/tree/main)

## Github File tree
```
├── Data
│   ├── metadata.csv
│   ├── Orthogroups.GeneCount.tsv
│   ├── Orthogroups.tsv
├── Data_Analysis_HTproject_files
│   ├── figure-markdown_strict
│   │   ├── unnamed-chunk-7-1.png
│   │   ├── unnamed-chunk-8-1.png
│   │   ├── unnamed-chunk-9-1.png
├── Data_Analysis_HTproject.md
├── Data_Analysis_HTproject.Rmd
├── Fusarium_HeatTolerance.Rproj
├── HPC_scripts
│   ├── Busco+basicstats.sh
│   ├── Clean_sort_mask.sh
│   ├── Orthofinder.sh
│   ├── predict.sh
├── README.md
├── Results
│   ├── core_orthogroups_all12.csv
│   ├── orthogroups_softcore_heat_tolerant.csv
│   ├── orthogroups_unique_heat_tolerant.csv
│   ├── orthogroups_unique_non_heat_tolerant.csv
```

