# Identify antibodies binding the same epitope.

<img align="right" src="https://github.com/Aurelien-Pelissier/Ab-binding/blob/main/img/binder.png" width=400>


Antibodies are protective proteins produced by B cells in response to the presence of foreign pathogens. Advances in Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) technologies have considerably increased the amount of repertoire data that is available for analysis and improved our understanding of the dynamics of B-cell repertoires in both individuals and populations. 

This repository address the specific task of **finding convergent evolution to antigen specificity accross antibodies from different repertoires**. These binders are typically from a different ancestor clones and are thus challenging to identify.

&nbsp;

## Launch the analysis

This repositroy combines the CDRsim [1], paratype [2] and Abligity [3] method, which are implemented in `Ab_binding.py`. Importantly, these rely on ANARCI and parapred which have quite specific requirements to run. They only runs on Linux, and parapred requires an old version of Tensofrlow incompatible with the newest version of python. Thus, we advice to first create a new Anaconda environmenet:

	- conda create -n parapred python=3.6
	- conda activate parapred
	- pip install -r parapred/requirements.txt

Then install ANARCI and Levenshtein with:

	- pip install python-Levenshtein
	- conda install -c bioconda anarci

Importantly, if you want to use parapred and Ab-Ligity, You need to get full antibody structures by first submitting your sequences to Ab structure inference pipeline such as [Repertoire Builder](https://sysimm.org/rep_builder/) [4].

Then, tu run an analysis, you can simply run an example analysis with `run_example.py`

&nbsp;


## Example: Self-epitopes in Rheumatoid Arthritis

Here we show how the method in this repository can identify antibodies binding to the known epitopes of Collagen type II, a well characterized antigen involved in Rheumatoid Arthrisis. Since more than 40 antibodies binding to CII were previously characterized~, we can identify the epitope reactivity of many of the sequences in the immune repertoires from different organs. 

<img src="https://github.com/Aurelien-Pelissier/Ab-binding/blob/main/img/RAmice.png" width=800>

&nbsp;


## References
[//]: <> (This may be the most platform independent comment)

[1] Pelissier, Aurelien, et al. "Convergent Evolution and B-Cell Recirculation in Germinal Centers in a Human Lymph Node." BioRxiv (2022).

[2] Richardson, Eve, et al. "A computational method for immune repertoire mining that identifies novel binders from different clonotypes, demonstrated by identifying anti-pertussis toxoid antibodies." mAbs. Vol. 13. No. 1. Taylor & Francis, 2021.

[3] Wong, Wing Ki, et al. "Ab-Ligity: identifying sequence-dissimilar antibodies that bind to the same epitope." MAbs. Vol. 13. No. 1. Taylor & Francis, 2021.

[4] Schritt, Dimitri, et al. "Repertoire Builder: high-throughput structural modeling of B and T cell receptors." Molecular Systems Design & Engineering 4.4 (2019): 761-768.
