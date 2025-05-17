import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from chemcluster import mol_to_base64_img
from rdkit import Chem

def test_mol_to_base64_img():
    smiles = "CCO"  # Ethanol
    mol = Chem.MolFromSmiles(smiles)

    result = mol_to_base64_img(mol)

    assert isinstance(result, str)
    assert result.startswith("data:image/png;base64,")
    assert len(result) > 100 