# Docking with AutoDock Vina or Smina

This script automates the process of docking ligands to a receptor using either **AutoDock Vina** or **Smina**. It also handles the preparation of ligands and receptors, generating the necessary PDBQT files required for docking.

## Requirements

Before using this script, ensure the following dependencies are installed:

1. **RDKit**: Required for ligand preparation, including 3D structure generation and force field optimization.
   - Install rdkit via pip
       ```bash
       pip install rdkit
       ```

2. **AutoDock Vina**: A molecular docking tool for drug discovery, required for docking.
   - Install Vina from [here]([http://vina.scripps.edu/download.html](https://vina.scripps.edu/downloads/))

3. **Smina**: An open-source fork of AutoDock Vina, used as an alternative docking tool.
   - Install Smina from [here]([https://github.com/glennon/smina](https://sourceforge.net/projects/smina/))

4. **ADFR Suite**: A suite of tools for molecular docking.
   - Install ADFR Suite from [here](https://ccsb.scripps.edu/adfr/downloads/)

5. **Open Babel**: A chemical toolbox for converting between different molecular structure file formats.
   - Install Open Babel from [here](https://openbabel.github.io/docs/Installation/install.html)

### Python Libraries
- **Subprocess**: For executing system commands.
- **RDKit**: For molecule manipulation and preparation.
- **Argparse**: For handling command-line arguments.

To install these Python dependencies, run:
```bash
pip install rdkit
```

### Additional tips
Ensure that the executables for **AutoDock Vina** or **Smina** are accessible from your command line. The default paths are expected to be in the system's `PATH` or can be specified directly when running the script.

## Usage

### Command-line Arguments
The script can be executed via the command line with the following arguments:

```bash
python3 docking_script.py --receptor <path_to_receptor.pdb> --ligand <path_to_ligand.sdf> --center <x_center y_center z_center> --size <x_size y_size z_size> [options]
```

#### Required Arguments:
- `--receptor`: Path to the receptor file (PDB format).
- `--ligand`: Path to the ligand file (PDB, MOL2, or SDF format) or a ligand SMILES string.
- `--center`: The center of the docking box in the format `[x, y, z]`.
- `--size`: The size of the docking box in the format `[x, y, z]`.

#### Optional Arguments:
- `--out`: Output name for the docking results (default: `docked`).
- `--exhaustiveness`: Exhaustiveness of the docking search (default: `8`).
- `--cpu`: Number of CPU cores to use (default: `1`).
- `--ff`: Force field for ligand preparation (default: `mmff`).
- `--receptor_prepared`: Whether the receptor is already prepared (pdbqt format).
- `--ligand_prepared`: Whether the ligand is already prepared (pdbqt format).
- `--output_directory`: Directory to save docking results (default: `docked_output`).
- `--docking_program`: Choose either `vina` or `smina` for docking (default: `vina`).
- `--vina_path`: Path to the Vina executable.
- `--smina_path`: Path to the Smina executable.
- `--prepare_receptor_path`: Path to the `prepare_receptor` script from the ADFR suite.

### Example Command

```bash
python docking_script.py --receptor receptor.pdb --ligand ligand.sdf --center 0 0 0 --size 20 20 20 --out docked_results --docking_program vina --vina_path /path/to/vina
```

This command will perform docking using **AutoDock Vina** and save the results in the `docked_results.pdbqt` file.

## Output

The output will include:
- Docked ligand structures in PDBQT format in a separate folder in the same working directory
- sdf and pdbqt files of prepared ligand in a separate folder in the same working directory (not created when --ligand_prepared is set to True)
- pdbqt file of prepared receptor in a separated folder in the same working directory (not created when --receptor_prepared is set to True)
  
