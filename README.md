# Docking with AutoDock Vina or Smina

This script automates the process of docking ligands to a receptor using either **AutoDock Vina** or **Smina**. It also handles the preparation of ligands and receptors, generating the necessary PDBQT files required for docking.

## Requirements

Before using this script, ensure the following dependencies are installed:

1. **RDKit**: Required for ligand preparation, including 3D structure generation and force field optimization.
   - Install RDKit via pip:
       ```bash
       pip install rdkit
       ```

2. **AutoDock Vina**: A molecular docking tool for drug discovery, required for docking.
   - Install Vina from [here](https://vina.scripps.edu/downloads/)

3. **Smina**: An open-source fork of AutoDock Vina, used as an alternative docking tool.
   - Install Smina from [here](https://sourceforge.net/projects/smina/)

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

The process of file preparation and docking can be run either through the command line or through a simple interactive interface.

### 1. Interactive Usage (Command-line Interactive Mode)
This method guides the user through the ligand and receptor preparation, docking setup, and execution interactively. You will be prompted to input various parameters such as receptor and ligand file paths, docking box dimensions, and other options.

To run the interactive mode, use the following command:
```bash
python3 docking_interactive.py
```

#### Prompts
You will be prompted to enter the following details:
- **Receptor File Path** (PDB format or PDBQT format).
- **Ligand File Path** (PDB, MOL2, SDF format, or SMILES string).
- **Docking Box Center** (x, y, z coordinates).
- **Docking Box Size** (x, y, z dimensions).
- **Output File Name** (optional, default: `docked`).
- **Exhaustiveness** (default: 8).
- **Number of CPUs** (default: 1).
- **Force Field** (for ligand preparation: `mmff` or `uff` or default None).
- **Receptor and Ligand Prepared Status** (whether they are already in pdbqt format).
- **Docking Program** (choose between `vina` or `smina`).
- **Executable Paths** for Vina and Smina (if not in PATH).
- **Prepare Receptor Path** (path to the `prepare_receptor` script).

### 2. Command-line Arguments (Non-Interactive Mode)
The script can be executed via the command line with the following arguments:

```bash
python3 docking.py --receptor <path_to_receptor.pdb> --ligand <path_to_ligand.sdf> --center <x_center y_center z_center> --size <x_size y_size z_size> [options]
```

#### Required Arguments:
- `--receptor`: Path to the receptor file (PDB format or pdbqt format).
- `--ligand`: Path to the ligand file (PDB, MOL2, or SDF format) or a ligand SMILES string.
- `--center`: The center of the docking box in the format `[x, y, z]`.
- `--size`: The size of the docking box in the format `[x, y, z]`.

#### Optional Arguments:
- `--out`: Output name for the docking results (default: `docked`).
- `--exhaustiveness`: Exhaustiveness of the docking search (default: `8`).
- `--cpu`: Number of CPU cores to use (default: `1`).
- `--ff`: Force field for ligand preparation (default: None).
- `--receptor_prepared`: Whether the receptor is already prepared (pdbqt format).
- `--ligand_prepared`: Whether the ligand is already prepared (pdbqt format).
- `--output_directory`: Directory to save docking results (default: `docked_output`).
- `--docking_program`: Choose either `vina` or `smina` for docking (default: `vina`).
- `--vina_path`: Path to the Vina executable.
- `--smina_path`: Path to the Smina executable.
- `--prepare_receptor_path`: Path to the `prepare_receptor` script from the ADFR suite.

## Output

The output will include:
- Docked ligand structures in PDBQT format in a separate folder in the same working directory.
- SDF and PDBQT files of the prepared ligand in a separate folder in the same working directory (not created when `--ligand_prepared` is set to True).
- PDBQT file of the prepared receptor in a separate folder in the same working directory (not created when `--receptor_prepared` is set to True).

## Troubleshooting

- **Error: "prepare_receptor not found"**  
  Ensure `prepare_receptor` from AutoDockTools is correctly installed and in PATH.
  
- **`obabel: command not found`**  
  Install Open Babel:
    ```bash
    sudo apt install openbabel
    ```

- **Vina or Smina executable not found**  
  Provide the full path to the executable when prompted or add them to PATH.

## Customization

- **Default paths and parameters** can be customized within the script by modifying default values in the `dock_interactive()` function.

