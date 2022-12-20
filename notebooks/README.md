# Source code

## Notebooks
- 1_chempert_preprocessing.ipynb. Regenerate ChemPert data and process it to generate the necessary vectors for each compound to run the analysis. Generate Figure 2.
- 2_correlation_chempert.ipynb. Conduct the correlation analysis. Generate Figure 3
- 2.1_DEG_analysis.ipynb. Run correlation analysis for section 3.1.1.
- 3_pathway_info_for_target_vectors.ipynb. Prepare the dataset to run PPI and pathway correlation analysis.
- 3.1_random_pathway_target_vectors.ipynb. Prepare the dataset to run pathway correlation analysis with randomly generated pathways.
- 3.2_random_network_target_vectors.ipynb. Prepare the dataset to run pathway correlation analysis with randomly generated networks.
- 4_pathway_correlation_analysis.ipynb. Run pathway / PPI correlation analysis.
- 4.1_random_networks_analysis.ipynb. Run correlation analysis for random networks. 
- 4.2_random_pathways_analysis.ipynb. Run correlation analysis for random pathways.
- 4.3_multifactorial_analysis.ipynb. Run correlation analysis on mjultifactorial pathways. 
- 5_similar_transcriptomic_profiles.ipynb. 3.3.2 analysis (Compounds targeting the same targets typically induce disparate transcriptomic responses). Generate Figure 4.
- 6_Network_figure.ipynb. Figure 5 of the manuscript.

## Scripts
- utils.py. Shared functions used in the analysis.
- viz.py. Code to generated Figure 5 of the manuscript.
