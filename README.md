# LassoPred

LassoPred predicts 3D structures using provided lasso peptide sequence data. It generates three possible optimized structures, prediction information, and MD input files.

Lasso peptide predictions for sequences within 30 residues take approximately 10 minutes.

## Requirements

- `joblib==1.4.2`
- `numpy==1.26.4`
- `scikit-learn==1.3.2`
- `scipy==1.13.0`
- `threadpoolctl==3.5.0`

Ensure your Linux environment (python=3.10) includes the following software:
- **Amber:** `pmemd.MPI`
- **PyMOL:** `pymol -c`

Hardware and parallel computing:

- **CPU:** A 16-core CPU is recommended.
- **MPI:** Install MPI and use the `mpirun` command to enable parallel computing capabilities.

## Install

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


## Usage

```bash
python script_name.py -seq <sequence> -tname <target_name> [-fdir <fold_direction>]
```

### Arguments
`-seq`, `--sequence`: The sequence for prediction (required).

`-tname`,` --target_name`: The target directory name (required).

`-fdir`, `--fold_direction`: Direction of lasso peptide ring folding (optional). Choices are left or right (default: right).

## Example

```bash
python predict.py -seq VGCETQEEVDELWAKLTADGGQEQPCAWLK -tname my_test1
python predict.py -seq GQIYDHPEVGIGAYGCEGLQR -tname my_test2
python predict.py -seq SIEDGTIKEAGSSQYYFV -tname my_test3
python predict.py -seq SVDDEGLREPMQPQLAWF -tname my_test4
python predict.py -seq SGIGDVFPEPNMVRRWD -tname my_test5
```

## Output

The directory under the target name you specified will include:

- Optimized structures: `min1.pdb`, `min2.pdb`, and `min3.pdb`. They are ordered by the predicted loop length.
- Prediction summary: `summary.csv`. This file records the predicted annotations.
- MD simulation files: `MDfiles1`, `MDfiles2`, `MDfiles3`. These folders include the force field files (`.prmtop` and `.inpcrd`) and correspond to the optimized structures.

## Note

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

