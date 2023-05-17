# Ab-Ligity: Identifying sequence-dissimilar antibodies that bind to the same epitope

Ab-Ligity describes the paratope or epitope similarity on antibody and antigen respectively. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine 
for development and testing purposes. See deployment for notes on how to deploy the project 
on a live system.

### Prerequisites

```
Biopython
Boost 1.66
GCC compiler
```
If you only have the antibody sequences / are working with antibody models:
```
ABodyBuilder or antibody modelling tools
Parapred: https://github.com/eliberis/parapred
```

### Installing

```
bash run.sh
```
This script compiles the C++ codes and run a test case.

## Running the tests

```
python3.6 abligity.py --mode 0 --pdb_dir example --pip_dir example -v
```

## Deployment

Add additional notes about how to deploy this on a live system

## Authors

* **Wing Ki Wong** - *Initial work*
* **Oxford Protein Informatics Group**

## License

This project is licensed under the BSD 3-clause License - see the [LICENSE.md](LICENSE.md) file for details.
You are using Parapred and ANARCI - their licenses are enclosed in the respective directories.

## Acknowledgments

* This work was supported by funding from the Engineering and Physical Sciences Research Council (EPSRC) 
and the Medical Research Council (MRC) [grant number EP/L016044/1].
