# Cell-Expression-Atlas-of-the-Autonomic-Nervous-System
Cell Expression Atlas of the Autonomic Nervous System

This app contains data from the following GEO accession ID's: 
# Sympathetic datasets
Superior cervical ganglion (GSE175421, GSE231766)
Stellate ganglion (GSE231924)
Celiac-superior mesenteric ganglion (GSE278457)

# Parasympathetic datasets
Dorsal motor nucleus of the vagus (GSE172411)
Nucleus ambiguus (GSE198709, GSE202760, GSE211538)

This dash app provides a platform for researchers to explore gene expression patterns in the mammalian autonomic nervous system (ANS). This “cell atlas” contains genetic data from over 70,000 cells across 8 different studies. It provides an integrated view of the autonomic nervous system facilitating future study on the cellular organization of the ANS.

The data is represented as 2 interactive 3-Dimentional dynamic UMAP scatter plots where each dot represents a single cell. Cells are colored by gene expression level on the left and dataset on the right. Expression data from over 20,000 genes can be searched using the search bar.

Various filters such as ANS division and ganglia are available for comparative analysis. In addition to the filters, there are two toggles: one for plotting neurons vs all cells, and one for switching between the original UMAP embedding and Harmony-integrated UMAP embedding. (Harmony integration reduces batch-specific biases and aligns datasets to reveal biologically relevant clusters.)

The Dash app is designed to be plug-and-play, so other labs can adapt the framework to work with their desired single-cell RNA-seq data (.h5ad format) and adata.obs metadata filters.

# Our Review
For more information on the ANS see our review paper at https://rdcu.be/euEhn
Code used for figure generation and raw data processing can be found on our paper-specific repository at https://github.com/YTwTJ/Molecular-and-Functional-Diversity-of-the-Autonomic-Nervous-System
