# Identifying antibodies binding the same epitope.

<img align="right" src="https://github.com/Aurelien-Pelissier/Ab-binding/blob/main/img/binder.png" width=400>


Antibodies are protective proteins produced by B-cells in response to the presence of foreign pathogens. Advances in Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) technologies have considerably increased the amount of repertoire data that is available for analysis and improved our understanding of the dynamics of B-cell repertoires in both individuals and populations. 

This repository address the specific task of **finding convergent specificity evolution to antigen accross antibodies from different repertoires**. These binders are typically from a different ancestor and are thus challenging to identify (different germline genes or CDR3 length).

### This repository support our publication:

[//]: <> (Pelissier, Aurelien, et al. "DOrnatien RHeumatoid arthirsis ACB model." BioRxiv 2022)

[[1]](https://www.life-science-alliance.org/content/6/11/e202301959) Pelissier, Aurelien, et al. "Convergent Evolution and B-Cell Recirculation in Germinal Centers in a Human Lymph Node." Life Science Alliance (2023).


&nbsp;


## Metrics

This repository contains script to combine different metrics comparing antibodies. It assess metrics based on functional similarities such as CDRsim [[1]](https://www.biorxiv.org/content/10.1101/2022.11.09.463832v9), Paratype [[2]](https://www.tandfonline.com/doi/full/10.1080/19420862.2020.1869406) and AbLigity [[3]](https://www.tandfonline.com/doi/full/10.1080/19420862.2021.1873478). Additionally, it also contains scripts to compute structure based metrics such as the tm-score and RMSD between antibody structures. 

## Prerequies

Depending on you metrics of interrest, you will need different requirements. First, most metrics rely on [parapred](https://github.com/eliberis/parapred) and [ANARCI](https://github.com/oxpig/ANARCI) which have quite specific requirements to run. They can only be executed on a **Linux** system, and parapred requires an old version of Tensorflow incompatible with the newest version of python. Thus, we advice to first create a new Anaconda Linux environment:

	- conda create -n parapred python=3.6
	- conda activate parapred
	- pip install -r tools/parapred/requirements.txt

Then installing ANARCI, Levenshtein and tmscoring with:

	- pip install python-Levenshtein
	- conda install -c bioconda anarci
	- pip install tmscoring



Furthermore, any structure based metrics, such as TM-score, RMSE or Ab-ligity, requires the pdb structure of all sequences you whish to compare. You can obtain these structures from your sequences with a homology modeling algorithm such as [Repertoire Builder](https://sysimm.org/rep_builder/) or recent with deep learning frameworks such as [IgFold](https://www.nature.com/articles/s41467-023-38063-x). 

## Launch the analysis

For a quick pairwise computation of the distance matrix of different metrics, you can simply run `python3 run_example.py`. 
However, to adapt your own sequences to the pipeline, you will need to run the `preprocess_data()` function in order to compute the CDRs and paratopes of your sequences (currently commented out in the file).
See the `raw_data/sample_repertoire.csv` for a file template to prepare your data for processing.

&nbsp;

## Do I need both chains ?

While the three methods are more reliable when paired heavy chain and light chain sequence information is available, they were shown to still perform well with the heavy chain information only. To get a higher precision, you can combine two of the methods together, but note that this will come at the cost of a lower recall. Check out [this study](https://www.biorxiv.org/content/biorxiv/early/2022/12/17/2022.11.09.463832/DC1/embed/media-1.pdf) (Section 5) for a benchmarking of the three methods on a Pertussis toxin binders dataset (Section 5 - Identifying common binders).

## When can we consider that two antibodies bind the same epitope ?
A similar paratope and CDR set of residue is a good indication that antibodies bind a similar epitope, even when not all CDRs have the same lenght. Independantly wether you used both chains or the heavy chain only, we advice the following thershold to determine if two antibodies bind the same epitope: 0.77 of CDRsim, 0.68 for paratype, and 0.18 for AbLigity (optimised from antibodies in the SabDab database). If your pair of antibodies is above the threshold for two of these metrics, there is a high chance that these antibodies share similar binding properties.


## Example: Self-epitopes in Rheumatoid Arthritis

Here we show how the method in this repository can identify antibodies binding to the known epitopes of Collagen type II, a well characterized antigen involved in Rheumatoid Arthrisis. Since more than 40 antibodies binding to CII on different epitopes were previously characterized (available in `known_binders.csv`), we can identify the epitope reactivity of many of the sequences in the immune repertoires, and observe the variability of the immune repertoire specificity within and accross organs. The figure bellow is taken from arthritis induced mice at different days after CII immunization.

<img src="https://github.com/Aurelien-Pelissier/Ab-binding/blob/main/img/RAmice.png" width=800>

&nbsp;


## References
[[2]](https://www.tandfonline.com/doi/full/10.1080/19420862.2020.1869406) Richardson, Eve, et al. "A computational method for immune repertoire mining that identifies novel binders from different clonotypes, demonstrated by identifying anti-pertussis toxoid antibodies." mAbs. Vol. 13. No. 1. Taylor & Francis, 2021.

[[3]](https://www.tandfonline.com/doi/full/10.1080/19420862.2021.1873478) Wong, Wing Ki, et al. "Ab-Ligity: identifying sequence-dissimilar antibodies that bind to the same epitope." MAbs. Vol. 13. No. 1. Taylor & Francis, 2021.
