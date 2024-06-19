# LassoPred

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for use](#instructions-for-use)
- [Note](#note)

# Overview

LassoPred predicts 3D structures using provided lasso peptide sequence data. It generates three possible optimized structures, prediction information, and MD input files.

Lasso peptide predictions for sequences within 30 residues take approximately 10 minutes.

We also provide a web server that can generate predicted structures including MD files. User can download the results and run the relax step on their server. [LassoPred_web](lassopred.accre.vanderbilt.edu)


# System Requirements

## Hardware Requirements
The `LassoPred` package requires a standard computer with sufficient CPU and GPU resources to support the operations defined by the user. For minimal performance, an 8-core CPU is recommended. For optimal performance, we recommend the following specifications::

- **CPU:** 8-core CPU (required for minimal performance)
- **GPU:** (optimal) NVIDIA GPU with at least 4GB VRAM, supporting CUDA Compute Capability 3.5 or higher
- **Memory:** 16GB of RAM

## Software requirements

### OS Requirements

The package development version is tested on *Linux* operating systems. The developmental version of the package has been tested on the following systems:

- Linux: CentOS 7.9 
- Mac OS

### Python Dependencies

- `joblib==1.4.2`
- `numpy==1.26.4`
- `scikit-learn==1.3.2`
- `scipy==1.13.0`
- `threadpoolctl==3.5.0`

### Additional Software Requirements

Ensure your Linux environment (python=3.10) includes the following software:
- Amber: `pmemd.MPI`
- PyMOL: `pymol -c`
- MPI: `mpirun`

# Installation Guide

## Setting up the environment

We recommend using a conda environment.

1. Build and activate a new environment:

```bash
conda create -n lassopred python=3.10
source activate lassopred
```

2. Install requirements:

```bash
pip install -r requirements.txt
```

3. Install PyMOL:

```bash
conda install conda-forge::pymol-open-source
```

4. Install Amber 

See [Amber Installation Guide](https://ambermd.org/Installation.php)

With the recommended specifications, the setup (except for Amber) will complete in about 10 minutes.

# Demo

Here we have several lasso peptide sequences that have not been reported in the [PDB](https://www.rcsb.org/). You can run LassoPred as follows:

```bash
python predict.py -seq VGCETQEEVDELWAKLTADGGQEQPCAWLK -tname my_test1
python predict.py -seq GQIYDHPEVGIGAYGCEGLQR -tname my_test2
python predict.py -seq SIEDGTIKEAGSSQYYFV -tname my_test3
python predict.py -seq SVDDEGLREPMQPQLAWF -tname my_test4
python predict.py -seq SGIGDVFPEPNMVRRWD -tname my_test5
```

The expected output files are included in the demos folder.

For an 8-core CPU computer, the average runtime for the samples above (13-30 amino acids) is approximately 4-12 minutes each.

# Instructions for use

LassoPred allows users to predict the 3D structure of a lasso peptide. We assume that the sequence has already been identified as a lasso peptide. If not, the sequence will be rejected, and errors may occur.

## Usage

To run the software on your data, use the following command:

```bash
python script_name.py -seq <sequence> -tname <target_name> [-fdir <fold_direction>]
```

## Arguments

`-seq`, `--sequence`: The sequence for prediction (required).

`-tname`,` --target_name`: The target directory name (required).

`-fdir`, `--fold_direction`: Direction of lasso peptide ring folding (optional). Choices are left or right (default: right).

## Reproduction Instructions

To reproduce all the quantitative results in the manuscript, follow these steps:

- Set up your environment as described in the "Setting up the environment" section.
- Use the sequences provided in the demo or your own sequences identified as lasso peptides.
- Run the predictions using the command mentioned in the "Usage" section.
- Compare the output files with the expected results included in the demos folder.

These instructions should enable you to replicate the results accurat

# Results

The directory under the target name you specified will include:

- Optimized structures: `min1.pdb`, `min2.pdb`, and `min3.pdb`. They are ordered by the predicted loop length.
- Prediction summary: `summary.csv`. This file records the predicted annotations.
- MD simulation files: `MDfiles1`, `MDfiles2`, `MDfiles3`. These folders include the force field files (`.prmtop` and `.inpcrd`) and correspond to the optimized structures.

# Note

For certain reasons, there may **not** be three predicted results.

- Sequence Too Short:

```bash
python predict.py -seq EIVVAGDETSGT -tname my_test
```

For this 12-amino acid sequence, with a minimum ring length of 7, minimum loop length of 3, and minimum tail length of 2, there will be only one predicted structure.

- Sequence Too Long:

```bash
python predict.py -seq VGCETQEEVDELWAKLTADGGQEQPCAWLKDKFGLSWQIVPRQLGELLSDPDPEKSQRVMQAMLQMSKIDIATLQAAYDGV -tname my_test
```

For this sequence, although our classifier predicts the top 3 loop lengths to be 49, 70, and 69, only one optimized structure and MD files will be generated for the first ranked prediction. This is because our scaffold supports loop lengths from 3 to 50.

- Sequence availability: 

**Valid Characters:** Ensure the sequence only contains standard amino acids (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y). If the sequence contains invalid characters, it will be rejected.

**Minimum Length:** The sequence must be at least 12 amino acids long. Sequences shorter than this will be rejected.

**Critical Positions:** The 7th, 8th, and 9th positions must contain at least one Glu (E) or Asp (D). If none of these positions contain E or D, the sequence will be rejected as no isopeptide bond can be predicted.

**Loop and Tail Lengths:** The loop length and tail length should be at least 5 residues combined. The shortest naturally discovered lasso peptide has a loop of 3 and a tail of 2. If the remaining length after the potential isopeptide position is less than 5, the sequence will be rejected.


## Citing this work

```bibtex
@article{ouyang2024predicting,
  title={Predicting 3D Structures of Lasso Peptides},
  author={Ouyang, Xingyu and Ran, Xinchun and Xu, Han and Zhao, Yi-Lei and Link, A James and Yang, Zhongyue},
  year={2024}
}
```
