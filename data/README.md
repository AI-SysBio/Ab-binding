# Germinal Center B-Cell Data

## Overview

This repository contains sequence-level and clone-level data from germinal center B-cell samples.

The dataset includes heavy-chain sequences, their clone assignments, and read-count information across **20 samples**, corresponding to **10 germinal centers (GCs)** with **2 replicates per GC**.

For each sequence, read abundance is provided in two forms:

* raw read counts across samples,
* normalized read counts, adjusted by the total number of reads in each sample.

Clone assignments are provided using two clonotyping approaches:

* **`HAC_0.16`**, the main clone definition used for downstream analyses;
* **`junction`**, an exact junction-matching definition provided as a reference.

Only heavy-chain information is included in this dataset. Light-chain sequences are not provided.

## Files

The repository contains clone assignments and associated sequence-level information under two clone definitions:

1. **`HAC_0.16`**: the recommended clone definition, based on hierarchical agglomerative clustering;
2. **`junction`**: a reference clone definition based on exact junction sequence matching.

For each clone definition, both sequence-level and clone-level information are provided.

## Sequence-level files

The sequence-level files contain one row per observed heavy-chain sequence. Each sequence is associated with its read-count profile across the 20 samples and with a clone assignment under the corresponding clonotyping strategy.

For the main sequence-level file, the column **`clone_label`** corresponds to the clone assignment obtained using the **`HAC_0.16`** clonotyping strategy.

### Main columns

| Column        | Description                                                                                                                                           |
| ------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `clone_label` | Clone ID assigned using the recommended **`HAC_0.16`** clonotyping strategy.                                                                          |
| `pop`         | A list containing the raw number of reads for the sequence in each of the 20 samples.                                                                 |
| `norm_pop`    | A list containing the normalized read abundance for the sequence in each sample, obtained by normalizing by the total number of reads in that sample. |

Additional sequence-level columns may include sequence identifiers, V gene, J gene, and junction sequence.

## Clone-level files

The clone-level files summarize information for each inferred clone. Each row corresponds to one clone and includes the clone ID together with associated clone-level metadata.

Unless otherwise specified, clone-level analyses should use the **`HAC_0.16`** clone definition.

## Clone definitions

### Hierarchical agglomerative clustering: `HAC_0.16`

This is the recommended clone definition.

Under this definition, two sequences are assigned to the same clone if they:

* use the same V gene,
* use the same J gene,
* and have at least 84% junction sequence identity.

This clone definition is used as the default for downstream analyses.

### Exact junction matching: `junction`

This reference definition assigns two sequences to the same clone only if they have the exact same junction sequence.

## Recommended usage

For downstream analyses, use the clone assignments from **`HAC_0.16`**. In the main sequence-level file, these assignments are stored in the **`clone_label`** column.

This clonotyping strategy provides the preferred balance between exact junction matching and more permissive clustering. The **`junction`** files are provided as a reference for comparison with exact junction-based clonotyping.

## Notes

* The dataset contains information from **20 samples** across **10 germinal centers**, with **2 replicates per GC**.
* The dataset contains **heavy-chain sequences only**; light-chain sequences are not included.
* Read counts are provided both as raw counts in `pop` and as sample-normalized values in `norm_pop`.
* Clone IDs may differ between the **`HAC_0.16`** and **`junction`** clone definitions.
* Unless otherwise specified, analyses should use the **`HAC_0.16`** clone definition.
