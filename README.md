# On the correspondence between the transcriptomic response of a compound and its effects on its targets

## Table of Contents

* [General Info](#general-info)
* [Dataset](#dataset)
* [Code structure](#code)
* [Citation](#citation)


## General Info
This repository contains code and data for our paper, "On the correspondence between the transcriptomic response of a
compound and its effects on its targets" (Hart *et al.*, 2022), which investigates the agreement between
transcriptomic profiles and target information for over 2,000 chemical compounds.

## Dataset
The data comes from [ChemPert](https://chempert.uni.lu/) a tools containing mappings between chemical perturbation and
transcriptional response for non-cancer cells.

The ChemPert dataset and other used and generated data is available at Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7164118.svg)](https://doi.org/10.5281/zenodo.7164118).

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
    |-- 3.1_random_pathway_target_vectors.ipynb
    |-- 3.2_random_network_target_vectors.ipynb
    |-- 4_pathway_correlation_analysis.ipynb
    |-- 4.1_random_networks_analysis.ipynb
    |-- 4.2_random_pathways_analysis.ipynb
    |-- 4.3_multifactorial_analysis.ipynb
    |-- 4.4_DEG_analysis.ipynb
    |-- 5_similar_transcriptomic_profiles.ipynb
    |-- 6_Network_figure.ipynb
    |-- utils.py
    `-- viz.py
```
