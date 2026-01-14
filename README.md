# mammalian_brain_consensus_taxonomy

Contains all code used for analysis in Werner, Suresh, French, and Gillis. The misplaced mouse Pax6 interneuron subclass: A cross-species transcriptomic reassignment. 2026.


Current analysis code can be found in code/jan_2025
Below is a breif description of each notebook. A detailed description of all processed data and the underlying data for each figure panel will be made avilable prior to publication.

orthology.Rmd – Identification of the 1to1-orthologs that are present across all the mouse and primate datasets

gaba_dataset_prep.Rmd – Initial processing of all mouse and primate datasets. Filtering for the CGE and MGE cells from the mouse whole brain atlases. 

gaba_metaneighbor_primate_consensus_mouse.Rmd – MetaNeighbor alignment of mouse clusters and primate subclasses. Includes the initial mouse MetaNeighbor assessments to identify the replicable mouse clusters. Also includes the mouse integration of the replicable clusters with high homology to primate subclasses. Plots for Figure 2 A, B, D, F, G and Supp. Fig. 4 A, B, C, and Supp. Fig. 6A

gaba_meta_markers.Rmd – Generation of mouse MetaMarkers using either original subclass labels or the homologous primate labels. Generation of the cross-primate MetaMarkers for just the consensus inhibitory neuron clusters. Generation of cross-mouse-primate MetaMarkers using the homologous cell-type annotations.

comparing_taxonomies.Rmd – Generation of cell-type taxonomies of the mouse datasets, comparing the original and homologous subclass labels. Plots for Figure 2 E and H, Supp. Fig. 8 A-I.

mouse_primate_subclass_marker_variability.Rmd – Assessments of primate subclass marker expression in mouse data using the original mouse subclass labels (namesake marker genes). Also includes the targeted assessments of marker enrichment and MetaNeighbor scores for the mouse Lamp5, Vip, and Sncg original clusters that map homologously to primate Lamp5, Vip, Pax6, and Sncg. Includes assessments of DE stats computed over increasingly diverged species. Includes the pareto fronts of markers across species. Plots for Figure 1C, E, F, and Figure 2C, J, K, L, and Figure 6A (mouse plots), and Supp. Fig. 7 A-E, and Supp. Fig. 9 A-C

allen_primate_metaMarker_exp.Rmd – UMAP and marker expression visualizations for the cross-primate Allen MTG data. Plots for Supp. Fig. 2A-D and Supp. Fig. 3A-D

synaptic_go_term_DE.Rmd - Identifying the GO terms with ion channels genes and the intersection with the mouse datasets. Integration of the mouse CGE cells using just ion channel genes and associated DE of those genes across the CGE subclasses. Plots for Fig. 3 A, B, and Supp. Fig. 10A, and Supp table 2

patch_seq_data.Rmd - Analysis of electrophysiological data from Patch-seq data. Adding homologous subclass labels to the Patch-seq data through our identified homologous labels in the Tasic 2018 dataset. Computation of omega-squared analysis and cell-cell similarity analysis for ephys and transcriptomic data. Identification of ephys features that distinguish the Pax6 cells. Plots for Figure 3 C-H and Supp. Fig. 10 B-D.

mouse_spatial_data.Rmd - Analysis for the top 100 nearest neighbors and layer enrichment for the mouse homologous subclass interneurons. Plots for Figure 4 A, B, C.

non_allen_primate_data.Rmd – Initial processing of the Yale and UT SouthWestern datasets (integration, de novo clustering) and quantification of primate and mouse subclass marker enrichments. Plots for Figure 1D, Figure 5 A-H, Figure 6A Yale and UT SouthWestern plots, and Supp. Fig. 5 A-D, G-J

non_allen_primate_metaNeighbor.Rmd – Primate and mouse MetaNeighbor assessments with the Yale and UT SouthWestern primate datasets. Plots for Supp. Fig. 5 E, F, K, L

tree_shrew_processing.Rmd – All initial processing of the tree-shrew dataset (filtering, identification of CGE and MGE cells, integration, de novo clustering, and subclass labeling) and all MetaNeighbor and marker enrichment assessments. Plots for Figure 6A tree-shrew plots, and Supp. Fig. 11 A-K

rat_processing.Rmd - All initial processing of the rat datasets, identification of CGE and MGE cells, de novo clustering, and subclass labeling. Plots for Supp. Fig. 12 A-J. 
