import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from chemcluster import show_3d_molecule
from rdkit import Chem
import py3Dmol

def test_show_3d_molecule():
    smiles = "CCO"  # Ethanol
    mol = Chem.MolFromSmiles(smiles)
    
    # Run the function
    viewer = show_3d_molecule(mol)
    
    assert viewer is not None
    assert isinstance(viewer, py3Dmol.view)