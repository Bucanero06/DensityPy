import os

import msgspec
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import AllChem


class DownloadSDFFileOutput(msgspec.Struct):
    path: str
    compound: pcp.Compound

def sdf_to_xyz(sdf_filename, output_folder):
    assert sdf_filename.endswith('.sdf'), f"Invalid SDF file: {sdf_filename}, must end with .sdf"
    # Verify/Create output folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Load SDF file
    suppl = Chem.SDMolSupplier(sdf_filename)
    print(f'{len(suppl)} molecules loaded from {sdf_filename}.')
    for idx, mol in enumerate(suppl):
        mol: Chem.Mol = mol
        # int_cid = mol.GetIntProp("PUBCHEM_COMPOUND_CID")
        # print(f'\n    {idx = }\n    {int_cid = }\n    {Chem.MolToSmiles(mol) = }\n')
        if mol is None:
            print(f"Skipping invalid molecule at index {idx} in {sdf_filename}.")
            continue

        try:
            # Add hydrogens
            mol_with_h = Chem.AddHs(mol)

            # Compute 3D coordinates for all atoms, including hydrogens
            AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol_with_h)

            # Extract molecule information
            atoms = mol_with_h.GetAtoms()
            num_atoms = len(atoms)

            # Generate XYZ content
            xyz_content = f"{num_atoms}\n"
            xyz_content += f"Molecule {idx} from {sdf_filename}\n"
            for atom in atoms:
                pos = mol_with_h.GetConformer().GetAtomPosition(atom.GetIdx())
                xyz_content += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"

            # Write to XYZ file
            name_of_file = sdf_filename.split('/')[-1].split('.')[0]
            xyz_filepath = os.path.join(output_folder, f"molecule_{idx}_from_{name_of_file}.xyz")
            with open(xyz_filepath, 'w') as xyz_file:
                xyz_file.write(xyz_content)
            print(f"Molecule {idx} processed successfully. Output saved to {xyz_filepath}")
        except Exception as e:
            print(f"Failed to process molecule {idx} due to an error: {e}")

def download_sdf_file(molecule_identifier, output_dir) -> DownloadSDFFileOutput:

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Function to download SDF
    def get_sdf(cid):
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d'

        response = requests.get(url)
        return response.text

    # Get PubChem data
    if molecule_identifier.isnumeric():
        compound = pcp.get_compounds(int(molecule_identifier), 'cid')[0]
    else:
        compound = pcp.get_compounds(molecule_identifier, 'name')[0]

    name = compound.iupac_name

    print(f"Downloading SDF for {molecule_identifier} (CID: {compound.cid})"
          f"\nName: {compound.iupac_name},\nFormula: {compound.molecular_formula},\nWeight: {compound.molecular_weight}")

    # Download structures
    sdf = get_sdf(compound.cid)
    output_path = f'{output_dir}/{name}.sdf'
    with open(output_path, 'w') as f:
        f.write(sdf)

    return DownloadSDFFileOutput(path=output_path, compound=compound)

if __name__ == '__main__':
    # Example usage
    molecules = '11465'
    output_directory = 'coordinate_files'
    downloaded_mol_info = download_sdf_file(molecules, output_directory)
    sdf_to_xyz(downloaded_mol_info.path, output_directory)
