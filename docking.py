#!/usr/bin/python3
import os
import subprocess
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def ligand_prep(lig,output_directory='ligand_output',ligand_name='ligand',ff=None):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    if lig.endswith('.sdf') or lig.endswith('.mol'):
        mol=Chem.MolFromMolFile(lig,removeHs=False)
    elif lig.endswith('.mol2'):
        mol=Chem.MolFromMol2File(lig)
    elif lig.endswith('.pdb'):
        mol=Chem.MolFromPDBFile(lig)
    elif lig.endswith('.smiles') or isinstance(lig,str):
        mol=Chem.MolFromSmiles(lig)
        mol=Chem.AddHs(mol)
    else:
        raise ValueError("Unsupported ligand format. Accepts .sdf, .mol2, .mol, .pdb, or SMILES format.")
    
    if mol is None:
        raise ValueError(f"Failed to read ligand from {lig}")
    
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol,AllChem.ETKDG())
    
    if ff is not None:
    
        if ff.lower() == 'mmff':
            AllChem.MMFFOptimizeMolecule(mol)
        elif ff.lower() == 'uff':
            AllChem.UFFOptimizeMolecule(mol)
        else:
            raise ValueError("Unsupported force field. Choose 'mmff' or 'uff'.")

    sdf_path=os.path.join(output_directory,f"{ligand_name}.sdf")
    with Chem.SDWriter(sdf_path) as writer:
        writer.write(mol)

    ligand_pdbqt=os.path.join(output_directory,f"{ligand_name}.pdbqt")
    obabel_cmd=f"obabel {sdf_path} -O {ligand_pdbqt} --gen3d --partialcharge"
    
    try:
        subprocess.run(obabel_cmd,shell=True,check=True)
    except subprocess.CalledProcessError as e:
        print("Error during ligand preparation.")
        print(e)

    return ligand_pdbqt


def receptor_prep(rec,output_directory='receptor_output',receptor_name='receptor',prepare_receptor_path='prepare_receptor'):
    if not rec.endswith('.pdb'):
        raise ValueError("Only pdb or pdbqt files are accepted")
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    receptor_output=os.path.join(output_directory,f'{receptor_name}.pdbqt')
    
    if not os.path.exists(prepare_receptor_path):
        raise ValueError(f"prepare_receptor not found at {prepare_receptor_path}. Please provide the correct path.")
    
    cleanup_command=f'{prepare_receptor_path} -r {rec} -o {receptor_output} -A hydrogens -U nphs_lps_waters_nonstdres -e'
    
    try:
        result=subprocess.run(cleanup_command,shell=True,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        result.check_returncode()
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error in receptor preparation. Check the input file or the 'prepare_receptor' tool. Error: {e.stderr.decode()}")
    
    return receptor_output


def dock(receptor,ligand,center,size,out='docked',exhaustiveness=8,cpu=1,ff=None,receptor_prepared=False,
         ligand_prepared=False,output_directory='docked_output',docking_program='vina',vina_path='vina', 
         smina_path='smina',prepare_receptor_path='prepare_receptor'):
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    if receptor_prepared:
        receptor_pdbqt=receptor
    else:
        receptor_pdbqt=receptor_prep(rec=receptor,prepare_receptor_path=prepare_receptor_path)
    
    if ligand_prepared:
        ligand_pdbqt=ligand
    else:
        ligand_pdbqt=ligand_prep(lig=ligand,ff=ff)
    
    docked_pdbqt=os.path.join(output_directory,f'{out}.pdbqt')
    docked_log=os.path.join(output_directory,f'{out}.log')
    
    if docking_program=='vina':
        if not os.path.exists(vina_path):
            raise ValueError(f"Vina executable not found at {vina_path}. Please provide the correct path.")
        
        command=f"{vina_path} --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} --out {docked_pdbqt} " \
                f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} " \
                f"--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} --exhaustiveness {exhaustiveness} --cpu {cpu}"
    
    elif docking_program=='smina':
        if not os.path.exists(smina_path):
            raise ValueError(f"Smina executable not found at {smina_path}. Please provide the correct path.")
        
        command=f"{smina_path} --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} --out {docked_pdbqt} " \
                f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} " \
                f"--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} --exhaustiveness {exhaustiveness} --cpu {cpu}"
    
    else:
        raise ValueError("Invalid docking program specified. Choose 'vina' or 'smina'.")
    
    try:
        subprocess.run(command,check=True,shell=True)
        print(f"Docking completed successfully. Results saved in {output_directory}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error during docking: {e.stderr.decode() if e.stderr else str(e)}")
    
    return docked_pdbqt


def main():
    parser=argparse.ArgumentParser(description="Docking with AutoDock Vina or Smina")
    parser.add_argument('--receptor',type=str,required=True,help="Receptor file (PDB format)")
    parser.add_argument('--ligand',type=str,required=True,help="Ligand file (PDB, MOL2, or SDF format) or Ligand Smiles String")
    parser.add_argument('--center',type=float,nargs=3,required=True,help="Center coordinates [x, y, z]")
    parser.add_argument('--size',type=float,nargs=3,required=True,help="Size of the docking box [x, y, z]")
    parser.add_argument('--out',type=str,default='docked',help="Output name for docking results")
    parser.add_argument('--exhaustiveness',type=int,default=8,help="Exhaustiveness of the search")
    parser.add_argument('--cpu',type=int,default=1,help="Number of CPU cores to use")
    parser.add_argument('--ff',type=str,default='mmff',help="Force field for ligand preparation")
    parser.add_argument('--receptor_prepared',type=bool,default=False,help="Whether receptor is already prepared (pdbqt format)")
    parser.add_argument('--ligand_prepared',type=bool,default=False,help="Whether ligand is already prepared (pdbqt format)")
    parser.add_argument('--output_directory',type=str,default='docked_output',help="Output directory for docking results")
    parser.add_argument('--docking_program',type=str,choices=['vina','smina'],default='vina',help="Docking program to use ('vina' or 'smina')")
    parser.add_argument('--vina_path',type=str,default='vina',help="Path to the vina executable")
    parser.add_argument('--smina_path',type=str,default='smina',help="Path to the smina executable")
    parser.add_argument('--prepare_receptor_path',type=str,default='prepare_receptor',help="Path to the prepare_receptor script from ADFR suite")

    args=parser.parse_args()

    dock(args.receptor,args.ligand,args.center,args.size,out=args.out,exhaustiveness=args.exhaustiveness,
         cpu=args.cpu,ff=args.ff,receptor_prepared=args.receptor_prepared,ligand_prepared=args.ligand_prepared,
         output_directory=args.output_directory,docking_program=args.docking_program,vina_path=args.vina_path,
         smina_path=args.smina_path,prepare_receptor_path=args.prepare_receptor_path)


if __name__ == "__main__":
    main()
