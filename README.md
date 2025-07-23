# mammalian_brain_consensus_taxonomy

Contains all code used for analysis in Werner, Suresh, Gillis. The misplaced mouse Pax6 interneuron subclass: A cross-species transcriptomic reassignment. 2025.


Current analysis code can be found in code/june_2025
Below is a breif description of each notebook. A detailed description of all processed data and the underlying data for each figure panel will be made avilable prior to publication.

gaba_dataset_prep.Rmd – Initial processing of all mouse and primate datasets

orthology.Rmd – Identification of the 1to1-orthologs that are present across all the mouse and primate datasets

gaba_metaneighbor_primate_consensus_mouse.Rmd – MetaNeighbor alignment of mouse clusters and primate subclasses. Includes the initial mouse MetaNeighbor assessments to identify the replicable mouse clusters. Also includes the mouse integration of the replicable clusters with high homology to primate subclasses.

gaba_meta_markers.Rmd – Generation of mouse MetaMarkers using either original subclass labels or the homologous primate labels. Generation of the cross-primate MetaMarkers for just the consensus inhibitory neuron clusters. Generation of cross-mouse-primate MetaMarkers using the homologous cell-type annotations.

mouse_primate_subclass_marker_variability.Rmd – Assessments of primate subclass marker expression in mouse data using the original mouse subclass labels. Also includes the targeted assessments of marker enrichment and MetaNeighbor scores for the mouse Lamp5, Vip, and Sncg original clusters that map homologously to primate Lamp5, Vip, Pax6, and Sncg.

comparing_taxonomies.Rmd – Generation of cell-type taxonomies of the mouse datasets, comparing the original and homologous subclass labels.

non_allen_primate_data.Rmd – Initial processing of the Yale and UT SouthWestern datasets (integration, de novo clustering) and quantification of primate and mouse subclass marker enrichments.

non_allen_primate_metaNeighbor.Rmd – Primate and mouse MetaNeighbor assessments with the Yale and UT SouthWestern primate datasets.

tree-shrew_processing.Rmd – All initial processing of the tree-shrew dataset (annotation, filtering, integration) and all MetaNeighbor and marker enrichment assessments.

allen_primate_metaMarker_exp.Rmd – UMAP and marker expression visualizations for the cross-primate Allen data.


