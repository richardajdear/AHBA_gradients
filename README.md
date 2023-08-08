# AHBA components
Code and data accompanying the article _Three components of human brain gene expression reflect normative developmental programmes with specific links to neurodevelopmental disorders_.

Python is used for analysis, R is used for plotting (ggplot+patchwork is just too good).

## Key outputs
Processing of the AHBA draws heavily on prior work by [Aurina Arnatkevičiūtė](https://www.sciencedirect.com/science/article/pii/S1053811919300114?via%3Dihub) and [Ross Markello](https://elifesciences.org/articles/72129).
- [processing_helpers.py](code/processing_helpers.py): this code applies the gene and region filters on top of the _abagen_ pipeline
- [ahba_dme_hcp_top8kgenes_scores.csv](outputs/ahba_dme_hcp_top8kgenes_scores.csv)
- [ahba_dme_hcp_top8kgenes_weights.csv](outputs/ahba_dme_hcp_top8kgenes_weights.csv)

## Analysis and figures
Figures and results are presented in jupyter notebooks that walk through the logic of each analysis. 
Figure-specific code is listed below each notebook, shared code is listed at the end.
- [fig1.ipynb](fig1.ipynb): generalisability summary & normative enrichments 
    - [triplets.py](code/triplets.py), [enrichments.py](code/enrichments.py), [enrichments_data.py](code/enrichments_data.py)
    - [fig1_plots.R](code/fig1_plots.R)
- [fig2.ipynb](fig2.ipynb): spatial associations to normative imaging
    - [maps_data.py](code/maps_data.py), [maps_analysis.py](code/maps_analysis.py), [maps_null_test.py](code/maps_null_test.py)
    - [fig2_plots.R](code/fig2_plots.R)
- [fig3.ipynb](fig3.ipynb): single-cell & developmental analyses
    - [single_cell.py](code/single_cell.py), [brainspan.py](code/brainspan.py)
    - [fig3_plots.R](code/fig3_plots.R)
- [fig4.ipynb](fig4.ipynb): disorder associations
    - [disorders_data.py](code/disorders_data.py), [disorders.py](code/disorders.py)
    - [fig4_plots.R](code/fig4_plots.R)
- [fig_extended.ipynb](fig_extended.ipynb): extended data & supplementary figures
    - [fig_extended.R](code/fig_extended.R)

Shared functions: these files are general depenencies for all the above
- [processing_helpers.py](code/processing_helpers.py): functions for processing the AHBA with abagen
- [gradientVersion.py](code/gradientVersion.py): class definition for PCA/DME object
- [analysis_helpers.py](code/analysis_helpers.py): general helper functions

Docker - to ease with installation of package dependencies do one of:
1. Pull the docker image richardajdear/ahba and run these analyses in a container (the docker image will automatically start a jupyter lab instance that you can access through a browser or an IDE like VScode)
2. Build your own docker image using the `Dockerfile`, or manually with the dependencies listed in `python-reqs.txt` and `r-reqs.R`