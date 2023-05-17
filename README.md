# Identify antibodies binding the same epitope.

<img align="right" src="https://github.com/Aurelien-Pelissier/Ab-binding/blob/main/img/binder.png" width=400>


Antibodies are protective proteins produced by B cells in response to the presence of foreign pathogens. Advances in Adaptive Immune Receptor Repertoire Sequencing (AIRR-Seq) technologies have considerably increased the amount of repertoire data that is available for analysis and improved our understanding of the dynamics of B-cell repertoires in both individuals and populations. 

**This repository address the specific task of finding convergent evolution to antigen specificity accross antibodies from different repertoires**. These binders are typically from a different ancestor clones and are thus challenging to identify.

&nbsp;

## Launch the analysis

This repositroy combines the CDRsim[1], paratype[2] and Abligity[3] method, which rely on ANARCI and parapred, which have specific requirements. Only runs on Linux, and parapred requires an old version of Tensofrlow so it's not compatible with the newest version of python.


## Example: Self-epitopes in Rheumatoid Arthritis

it can provide additional insights into the analysis of immune repertoires, such as finding convergent
epitope reactivity between independent samples. it can provide additional insights into the analysis of immune repertoires, such as finding convergent
epitope reactivity between independent samples. it can provide additional insights into the analysis of immune repertoires, such as finding convergent
epitope reactivity between independent samples. it can provide additional insights into the analysis of immune repertoires, such as finding convergent
epitope reactivity between independent samples. it can provide additional insights into the analysis of immune repertoires, such as finding convergent
epitope reactivity between independent samples. it can provide additional insights into the analysis of immune repertoires, such as finding convergent
epitope reactivity between independent samples. 

&nbsp;

<img src="https://github.com/Aurelien-Pelissier/Ab-binding/blob/main/img/RAmice.png" width=800>



## References
[//]: <> (This may be the most platform independent comment)

[1] Pelissier, Aurelien, et al. "Convergent Evolution and B-Cell Recirculation in Germinal Centers in a Human Lymph Node." BioRxiv (2022).

[2] Richardson, Eve, et al. "A computational method for immune repertoire mining that identifies novel binders from different clonotypes, demonstrated by identifying anti-pertussis toxoid antibodies." mAbs. Vol. 13. No. 1. Taylor & Francis, 2021.

[3] Wong, Wing Ki, et al. "Ab-Ligity: identifying sequence-dissimilar antibodies that bind to the same epitope." MAbs. Vol. 13. No. 1. Taylor & Francis, 2021.
