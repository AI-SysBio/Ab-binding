# Identifying antibodies binding the same epitope.

<img align="right" src="https://github.com/Aurelien-Pelissier/Ab-binding/blob/main/img/binder.png" width=400>


Antibodies are protective proteins produced by B-cells in response to the presence of foreign pathogens. Advances in Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) technologies have considerably increased the amount of repertoire data that is available for analysis and improved our understanding of the dynamics of B-cell repertoires in both individuals and populations. 

This repository address the specific task of **finding convergent specificity evolution to antigen accross antibodies from different repertoires**. These binders are typically from a different ancestor and are thus challenging to identify (different germline genes or CDR3 length).

&nbsp;

## Launch the analysis

The script in `Ab_binding.py` combines the CDRsim [[1]](https://www.biorxiv.org/content/10.1101/2022.11.09.463832v9), Paratype [[2]](https://www.tandfonline.com/doi/full/10.1080/19420862.2020.1869406) and AbLigity [[3]](https://www.tandfonline.com/doi/full/10.1080/19420862.2021.1873478) framework to identify potential common binders. Importantly, these rely on [ANARCI](https://github.com/oxpig/ANARCI) and [parapred](https://github.com/eliberis/parapred) which have quite specific requirements to run. They can only be executed on a **Linux** system, and parapred requires an old version of Tensorflow incompatible with the newest version of python. Thus, we advice to first create a new Anaconda Linux environment:

	- conda create -n parapred python=3.6
	- conda activate parapred
	- pip install -r parapred/requirements.txt

Then install ANARCI and Levenshtein with:

	- pip install python-Levenshtein
	- conda install -c bioconda anarci

Importantly, if you want to use AbLigity, You need to get full antibody structures by first submitting your sequences to Ab structure inference pipeline such as [Repertoire Builder](https://sysimm.org/rep_builder/). Then, to run an analysis, you can simply run an example analysis with `run_example.py`

&nbsp;

## Do I need both chains ?

While the three methods are more reliable when paired heavy chain and light chain sequence informationis available, they were shown to still perform well with the heavy chain information only.


## Example: Self-epitopes in Rheumatoid Arthritis

Here we show how the method in this repository can identify antibodies binding to the known epitopes of Collagen type II, a well characterized antigen involved in Rheumatoid Arthrisis. Since more than 40 antibodies binding to CII were previously characterized, we can identify the epitope reactivity of many of the sequences in the immune repertoires from different organs, and observe the variability of the immune repertoire specificity within and accross these organs. The figure bellow is taken from arthritis induced mice at different days after Collagen II immunization.

<img src="https://github.com/Aurelien-Pelissier/Ab-binding/blob/main/img/RAmice.png" width=800>

&nbsp;


## References
[//]: <> (This may be the most platform independent comment)

[[1]](https://www.biorxiv.org/content/10.1101/2022.11.09.463832v9) Pelissier, Aurelien, et al. "Convergent Evolution and B-Cell Recirculation in Germinal Centers in a Human Lymph Node." BioRxiv (2022).

[[2]](https://www.tandfonline.com/doi/full/10.1080/19420862.2020.1869406) Richardson, Eve, et al. "A computational method for immune repertoire mining that identifies novel binders from different clonotypes, demonstrated by identifying anti-pertussis toxoid antibodies." mAbs. Vol. 13. No. 1. Taylor & Francis, 2021.

[[3]](https://www.tandfonline.com/doi/full/10.1080/19420862.2021.1873478) Wong, Wing Ki, et al. "Ab-Ligity: identifying sequence-dissimilar antibodies that bind to the same epitope." MAbs. Vol. 13. No. 1. Taylor & Francis, 2021.
