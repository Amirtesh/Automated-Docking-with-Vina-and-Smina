#!/usr/bin/python3
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem


def ligand_prep(lig, output_directory='ligand_output', ligand_name='ligand', ff='mmff'):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    if lig.endswith('.sdf') or lig.endswith('.mol'):
        mol=Chem.MolFromMolFile(lig, removeHs=False)
    elif lig.endswith('.mol2'):
        mol=Chem.MolFromMol2File(lig)
    elif lig.endswith('.pdb'):
        mol=Chem.MolFromPDBFile(lig)
    elif lig.endswith('.smiles') or isinstance(lig, str):
        mol=Chem.MolFromSmiles(lig)
        mol=Chem.AddHs(mol)
    else:
        raise ValueError("Unsupported ligand format. Accepts .sdf, .mol2, .mol, .pdb, or SMILES format.")

    if mol is None:
        raise ValueError(f"Failed to read ligand from {lig}")

    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    if ff.lower() == 'mmff':
        AllChem.MMFFOptimizeMolecule(mol)
    elif ff.lower() == 'uff':
        AllChem.UFFOptimizeMolecule(mol)
    else:
        raise ValueError("Unsupported force field. Choose 'mmff' or 'uff'.")

    sdf_path=os.path.join(output_directory, f"{ligand_name}.sdf")
    with Chem.SDWriter(sdf_path) as writer:
        writer.write(mol)

    ligand_pdbqt=os.path.join(output_directory, f"{ligand_name}.pdbqt")
    obabel_cmd=f"obabel {sdf_path} -O {ligand_pdbqt} --gen3d --partialcharge"

    try:
        subprocess.run(obabel_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print("Error during ligand preparation.")
        print(e)

    return ligand_pdbqt



def receptor_prep(rec, output_directory='receptor_output', receptor_name='receptor', prepare_receptor_path='prepare_receptor'):
    if not rec.endswith('.pdb'):
        raise ValueError("Only pdb or pdbqt files are accepted")

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    receptor_output=os.path.join(output_directory, f'{receptor_name}.pdbqt')

    if not os.path.exists(prepare_receptor_path):
        raise ValueError(f"prepare_receptor not found at {prepare_receptor_path}. Please provide the correct path.")

    cleanup_command=f'{prepare_receptor_path} -r {rec} -o {receptor_output} -A hydrogens -U nphs_lps_waters_nonstdres -e'

    try:
        subprocess.run(cleanup_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error in receptor preparation: {e.stderr.decode()}")

    return receptor_output


def dock_interactive():
    print("\nWelcome to Interactive Docking Setup!\n")
    
    receptor=input("Enter receptor file path (PDB format or PDBQT format): ")
    ligand=input("Enter ligand file path (PDB, MOL2, SDF) or SMILES string: ")
    center_x=float(input("Enter docking box center (x): "))
    center_y=float(input("Enter docking box center (y): "))
    center_z=float(input("Enter docking box center (z): "))
    size_x=float(input("Enter docking box size (x): "))
    size_y=float(input("Enter docking box size (y): "))
    size_z=float(input("Enter docking box size (z): "))
    
    out=input("Enter output file name (default 'docked'): ") or 'docked'
    exhaustiveness=int(input("Enter exhaustiveness (default 8): ") or 8)
    cpu=int(input("Enter number of CPUs (default 1): ") or 1)
    ff=input("Force field for ligand preparation (mmff or uff, default mmff): ") or 'mmff'
    
    receptor_prepared=input("Is the receptor already prepared (pdbqt)? (y/n): ").strip().lower() == 'y'
    ligand_prepared=input("Is the ligand already prepared (pdbqt)? (y/n): ").strip().lower() == 'y'
    
    docking_program=input("Docking program (vina or smina, default vina): ") or 'vina'
    
    if docking_program == 'vina':
        vina_path=input("Path to Vina executable (default 'vina'): ") or 'vina'
        smina_path=None
    else:
        smina_path=input("Path to Smina executable (default 'smina'): ") or 'smina'
        vina_path=None
    
    prepare_receptor_path=input("Path to prepare_receptor (default 'prepare_receptor'): ") or 'prepare_receptor'
    output_directory=input("Output directory (default 'docked_output'): ") or 'docked_output'

    dock(
        receptor=receptor,
        ligand=ligand,
        center=[center_x, center_y, center_z],
        size=[size_x, size_y, size_z],
        out=out,
        exhaustiveness=exhaustiveness,
        cpu=cpu,
        ff=ff,
        receptor_prepared=receptor_prepared,
        ligand_prepared=ligand_prepared,
        docking_program=docking_program.lower().strip(),
        vina_path=vina_path,
        smina_path=smina_path,
        prepare_receptor_path=prepare_receptor_path,
        output_directory=output_directory
    )



def dock(receptor, ligand, center, size, out='docked', exhaustiveness=8, cpu=1, ff='mmff', receptor_prepared=False,
         ligand_prepared=False, output_directory='docked_output', docking_program='vina', vina_path='vina',
         smina_path='smina', prepare_receptor_path='prepare_receptor'):
    print("Starting docking...")

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    if receptor_prepared:
        receptor_pdbqt=receptor
    else:
        receptor_pdbqt=receptor_prep(receptor, prepare_receptor_path=prepare_receptor_path)

    if ligand_prepared:
        ligand_pdbqt=ligand
    else:
        ligand_pdbqt=ligand_prep(ligand, ff=ff)

    docked_pdbqt=os.path.join(output_directory, f'{out}.pdbqt')

    if docking_program=='vina':
        command=f"{vina_path} --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} --out {docked_pdbqt} " \
              f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} " \
              f"--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} --exhaustiveness {exhaustiveness} --cpu {cpu}"
    elif docking_program=='smina':
        command=f"{smina_path} --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} --out {docked_pdbqt} " \
              f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} " \
              f"--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} --exhaustiveness {exhaustiveness} --cpu {cpu}"
    else:
        raise ValueError("Only 'vina' or 'smina' is supported")

    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Docking completed successfully. Results saved in {output_directory}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error during docking: {e.stderr.decode() if e.stderr else str(e)}")


if __name__ == "__main__":
    dock_interactive()
