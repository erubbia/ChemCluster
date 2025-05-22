from rdkit import Chem
from rdkit.Chem import (
    Descriptors, Crippen, Lipinski, rdMolDescriptors,
    Draw, AllChem
)
import py3Dmol
from io import BytesIO
import base64


def calculate_properties(mol, mol_name="Unknown"):
    """ Calculates several molecular properties from an RDKit Mol object """
    return {
        "Molecule": mol_name,
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Crippen.MolLogP(mol), 2),
        "H-Bond Donors": Lipinski.NumHDonors(mol),
        "H-Bond Acceptors": Lipinski.NumHAcceptors(mol),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Aromatic Rings": Lipinski.NumAromaticRings(mol),
        "Heavy Atom Count": mol.GetNumHeavyAtoms()
    }

def get_fingerprint(mol):
    """ Generates a Morgan (circular) fingerprint bit vector for a given molecule """
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

def clean_smiles_list(smiles_list):
    """ Filters and converts a list of SMILES strings into RDKit Mol objects """
    mols, valid_smiles = [], []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols.append(mol)
            valid_smiles.append(smi)
    return mols, valid_smiles

def show_3d_molecule(mol, confId=-1):
    """  Displays an interactive 3D view of an RDKit molecule using py3Dmol """
    if not mol.GetNumConformers():
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mb = Chem.MolToMolBlock(mol, confId=confId)
    viewer = py3Dmol.view(width=300, height=300)
    viewer.addModel(mb, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.setBackgroundColor('0xffffff')
    viewer.zoomTo()
    return viewer

def mol_to_base64_img(mol):
    """ Converts an RDKit molecule to a base64-encoded PNG image """
    img = Draw.MolToImage(mol, size=(200, 200))
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    buffer.seek(0)
    return "data:image/png;base64," + base64.b64encode(buffer.read()).decode()
