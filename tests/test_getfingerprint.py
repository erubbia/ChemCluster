import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from chemcluster import get_fingerprint
from rdkit import Chem
from rdkit.DataStructs import ExplicitBitVect

def test_get_fingerprint():
    smiles = "CCO"  # ethanol
    mol = Chem.MolFromSmiles(smiles)
    fp = get_fingerprint(mol)

    assert isinstance(fp, ExplicitBitVect)
    assert fp.GetNumBits() == 2048
    assert fp.GetOnBits()  # There should be at least some bits turned on
    