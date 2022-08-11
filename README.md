# On the correspondence between the transcriptomic response of a compound and its effects on its targets

## Table of Contents

* [General Info](#general-info)
* [Dataset](#dataset)
* [Code structure](#code)
* [Citation](#citation)


## General Info
This repository contains code and data for our paper, "On the correspondence between the transcriptomic response of a
compound and its effects on its targets" (Engler *et al.*, 2022), which investigates the agreement between
transcriptomic profiles and target information for over 2,000 chemical compounds.

## Dataset
The data comes from [ChemPert](https://chempert.uni.lu/) a tools containing mappings between chemical perturbation and
transcriptional response for non-cancer cells.

## Code
The project has the following structure:

```
|-- data (main data folder)
|   |-- ChemPert_data
|   |-- Protein_classes
|   |-- Protein_mappings
|   |-- Protein_pathways
|   |-- Transcriptional_data_frames
|   |-- target_data_frames
|-- figures (paper figures)
`-- notebooks (correlation analyses)
    |-- 1_chempert_preprocessing.ipynb
    |-- 2_correlation_chempert.ipynb
    |-- 3_pathway_info_for_target_vectors.ipynb
    |-- 4_pathway_correlation_analysis.ipynb
    |-- 5_similar_transcriptomic_profiles.ipynb
    |-- 6_Network_figure.ipynb
    |-- utils.py
    `-- viz.py
```
