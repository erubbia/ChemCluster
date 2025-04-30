import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
import streamlit.components.v1 as components
from sklearn.decomposition import PCA
import plotly.express as px
from rdkit.Chem import AllChem, DataStructs

# Helper to calculate properties
def calculate_properties(mol, mol_name="Unknown"):
    properties = {
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
    return properties

# Helper to get fingerprint
def get_fingerprint(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

# Helper to calculate similarity matrix
def calculate_similarity_matrix(fingerprints):
    n = len(fingerprints)
    sim_matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                sim_matrix[i][j] = sim
    return sim_matrix

# Helper to plot PCA
def plot_pca(similarity_matrix, names):
    import numpy as np

    # Convert similarity to distance
    distance_matrix = 1 - np.array(similarity_matrix)

    pca = PCA(n_components=2)
    components_2d = pca.fit_transform(distance_matrix)

    plot_df = pd.DataFrame({
        "PCA 1": components_2d[:, 0],
        "PCA 2": components_2d[:, 1],
        "Molecule": names
    })

    fig = px.scatter(
        plot_df,
        x="PCA 1",
        y="PCA 2",
        hover_data=["Molecule"],
        color_discrete_sequence=["#4dabf7"] * len(plot_df),
        height=600,
        width=800
    )

    fig.update_traces(marker=dict(size=12, line=dict(width=1, color='DarkSlateGrey')))
    fig.update_layout(
        title="Molecule Similarity Map (Fingerprint PCA)",
        title_x=0.5,
        plot_bgcolor="#f7f9fb",
        paper_bgcolor="#f7f9fb",
        font=dict(size=14, color="#1f1f1f"),
        xaxis=dict(title="PCA Component 1", gridcolor="#d0d6e2", linecolor="#bcc2cd", ticks="outside"),
        yaxis=dict(title="PCA Component 2", gridcolor="#d0d6e2", linecolor="#bcc2cd", ticks="outside")
    )

    return fig

# App title
st.title("MolSearch: Explore Molecules")

st.write("""
Welcome to MolSearch! 
Upload molecules, input SMILES, draw structures, and explore chemical spaces interactively!
""")

input_method = st.radio(
    "How would you like to input your molecules:",
    ("Upload a file", "Input SMILES text", "Draw molecule with bonds (ChemDraw style)")
)

# File upload or SMILES input
mols = []
names = []
properties_list = []

if input_method == "Upload a file":
    uploaded_file = st.file_uploader("Upload your molecule file (.sdf, .mol, .csv with SMILES)", type=["sdf", "mol", "csv"])
    if uploaded_file is not None:
        st.success("File uploaded successfully! 🎉")
        file_type = uploaded_file.name.split(".")[-1].lower()

        if file_type == "sdf":
            suppl = Chem.ForwardSDMolSupplier(uploaded_file)
            mols = [mol for mol in suppl if mol is not None]

        elif file_type == "mol":
            mol = Chem.MolFromMolBlock(uploaded_file.read().decode("utf-8"))
            if mol:
                mols = [mol]

        elif file_type == "csv":
            df = pd.read_csv(uploaded_file)
            smiles_col = "SMILES" if "SMILES" in df.columns else "smiles"
            mols = [Chem.MolFromSmiles(sm) for sm in df[smiles_col].dropna()]
            names = df[smiles_col].tolist()

elif input_method == "Input SMILES text":
    smiles_input = st.text_area("Enter one or more SMILES (one per line):", height=200)
    if smiles_input:
        smiles_list = smiles_input.strip().splitlines()
        mols = [Chem.MolFromSmiles(smi) for smi in smiles_list if Chem.MolFromSmiles(smi)]
        names = smiles_list

elif input_method == "Draw molecule with bonds (ChemDraw style)":
    st.write("Draw your molecule below:")
    components.html(
        """
        <iframe src="https://partridgejiang.github.io/Kekule.js/demos/items/chemEditor/chemEditor.html" width="900" height="600" style="border:none;"></iframe>
        """,
        height=650,
    )
    st.info("After drawing, click 'Save/Export' inside the editor to get the SMILES, then paste it into the 'Input SMILES' section.")

# If molecules were loaded
if mols:
    fingerprints = [get_fingerprint(mol) for mol in mols]
    similarity_matrix = calculate_similarity_matrix(fingerprints)

    for i, mol in enumerate(mols):
        props = calculate_properties(mol, mol_name=names[i] if names else f"Molecule {i+1}")
        properties_list.append(props)

        col1, col2 = st.columns([1, 2])
        with col1:
            st.image(Draw.MolToImage(mol, size=(200, 200)), caption=f"Molecule {i+1}")
        with col2:
            st.write("**Properties:**")
            st.json(props)

            similarities = similarity_matrix[i]
            best_idx = similarities.index(max(similarities))
            best_score = similarities[best_idx]
            st.write(f"🔵 Most similar to: Molecule {best_idx+1} (Similarity {best_score:.2f})")

    df_props = pd.DataFrame(properties_list)
    st.write("📋 Preview of molecule properties:")
    st.dataframe(df_props)

    csv = df_props.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Download Properties as CSV", data=csv, file_name="molecule_properties.csv", mime="text/csv")

    if st.button("🔵 Show Similarity Map (Fingerprint PCA)"):
        fig = plot_pca(similarity_matrix, [prop['Molecule'] for prop in properties_list])
        st.plotly_chart(fig)
