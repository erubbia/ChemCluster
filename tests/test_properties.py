from rdkit import Chem
from chemcluster.app import calculate_properties  # adjust this import if needed

def test_calculate_properties():
    smiles = "CCO"  # ethanol
    mol = Chem.MolFromSmiles(smiles)
    result = calculate_properties(mol, mol_name="ethanol")

    assert isinstance(result, dict)
    assert result["Molecule"] == "ethanol"
    assert "Molecular Weight" in result
    assert "LogP" in result
    assert "H-Bond Donors" in result
    assert "H-Bond Acceptors" in result
    assert "TPSA" in result
    assert "Rotatable Bonds" in result
    assert "Aromatic Rings" in result
    assert "Heavy Atom Count" in result

    # Check known property value (approximate)
    assert abs(result["Molecular Weight"] - 46.07) < 0.1