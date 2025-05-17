
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from chemcluster import clean_smiles_list
from rdkit import Chem

def test_clean_smiles_list():
    smiles_input = ["CCO", "INVALID", "C1CCCCC1"]  # ethanol, invalid, cyclohexane
    mols, valid_smiles = clean_smiles_list(smiles_input)

    assert len(mols) == 2
    assert len(valid_smiles) == 2
    assert "CCO" in valid_smiles
    assert "C1CCCCC1" in valid_smiles
    assert all(isinstance(m, Chem.Mol) for m in mols)