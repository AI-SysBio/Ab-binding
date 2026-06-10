# Germinal Center B-Cell Data

## Overview

This repository contains sequence-level and clone-level data from germinal center B-cell samples.

The dataset includes all observed sequences, their clone assignments, and read-count information across **20 samples**, corresponding to **10 germinal centers (GCs)** with **2 replicates per GC**.

For each sequence, read abundance is provided in two forms:

* raw read counts across samples,
* normalized read counts, adjusted by the total number of reads in each sample.

Clones were defined using three alternative clonotyping strategies:

* exact junction sequence matching,
* hierarchical agglomerative clustering with 84% junction sequence identity, denoted **`HAC_0.16`**,
* hierarchical agglomerative clustering with 80% junction sequence identity, denoted **`HAC_0.2`**.

For downstream analyses, we recommend using **`HAC_0.16`** as the default clone definition.

## Files

The repository contains two main data files:

1. a **sequence-level file**, with one row per observed sequence;
2. a **clone-level file**, with one row per inferred clone.

## Sequence-level file

The sequence-level file contains one row per sequence. Each sequence is associated with its read-count profile across the 20 samples and with clone assignments under the three clonotyping strategies.

### Main columns

| Column     | Description                                                                                                                                                                                |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `pop`      | A list containing the raw number of reads for the sequence in each of the 20 samples.                                                                                                      |
| `norm_pop` | A list containing the normalized read abundance for the sequence in each sample, obtained by normalizing by the total number of reads in that sample.                                      |
| `junction` | Clone assignment based on exact junction sequence matching. Sequences are assigned to the same clone if they have identical junction sequences.                                            |
| `HAC_0.16` | Clone assignment based on hierarchical agglomerative clustering using the same V gene, same J gene, and at least 84% junction sequence identity. This is the recommended clone definition. |
| `HAC_0.2`  | Clone assignment based on hierarchical agglomerative clustering using the same V gene, same J gene, and at least 80% junction sequence identity.                                           |

## Clone-level file

The clone-level file summarizes information for each inferred clone. It contains clone IDs and associated clone-level metadata.

Unless otherwise specified, clone-level analyses should use the clone IDs defined by **`HAC_0.16`**.

## Clone definitions

Three clonotyping strategies were used.

### Exact junction matching: `junction`

Under this definition, two sequences are assigned to the same clone if they have the exact same junction sequence.

### Hierarchical agglomerative clustering: `HAC_0.16`

Under this definition, two sequences are assigned to the same clone if they:

* use the same V gene,
* use the same J gene,
* and have at least 84% junction sequence identity.

This is the recommended clone definition and should be used as the default for most analyses.

### Hierarchical agglomerative clustering: `HAC_0.2`

Under this definition, two sequences are assigned to the same clone if they:

* use the same V gene,
* use the same J gene,
* and have at least 80% junction sequence identity.

## Recommended usage

For downstream analyses, use the clone assignments in the **`HAC_0.16`** column. This clonotyping strategy provides the preferred balance between exact junction matching and more permissive clustering.

## Notes

* The dataset contains information from **20 samples** across **10 germinal centers**, with **2 replicates per GC**.
* Read counts are provided both as raw counts in `pop` and as sample-normalized values in `norm_pop`.
* Clone IDs may differ depending on the clonotyping strategy used.
* Unless otherwise specified, analyses should use the **`HAC_0.16`** clone definition.
