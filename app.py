
# --- MolSearch: Final Split Version (Single vs Dataset Mode) ---

import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, Draw, AllChem, DataStructs, rdDistGeom, rdMolAlign
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np
import plotly.express as px
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.ML.Cluster import Butina

# Aesthetic setup
st.set_page_config(page_title="MolSearch", layout="wide")
st.markdown("""
    <style>
    html, body, [class*="css"] {
        font-family: 'Nunito', sans-serif;
        font-size: 18px;
        background-color: #f9f9f9;
    }
    h1, h2, h3 {
        color: #1f3b57;
    }
    .stButton > button {
        background-color: #BCD4CC;
        color: #1f1f1f;
        border-radius: 10px;
        padding: 0.6em 1.2em;
        font-size: 18px;
    }
    .stButton > button:hover {
        background-color: #E1ADAD;
        color: #1f1f1f;
        transform: scale(1.05);
    }
    </style>
""", unsafe_allow_html=True)


def calculate_properties(mol, mol_name="Unknown"):
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
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

def clean_smiles_list(smiles_list):
    mols, valid_smiles = [], []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols.append(mol)
            valid_smiles.append(smi)
    return mols, valid_smiles

def show_3d_molecule(mol, confId=-1):
    if not mol.GetNumConformers():
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mb = Chem.MolToMolBlock(mol, confId=confId)
    viewer = py3Dmol.view(width=300, height=300)
    viewer.addModel(mb, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.setBackgroundColor('0xffffff')
    viewer.zoomTo()
    return viewer

# Title
st.markdown("""
<h1><span style='color:#1f3b57;'>Mol</span><span style='color:#BCD4CC;'>Search</span></h1>
""", unsafe_allow_html=True)

mode = st.radio("Choose analysis mode:", ["📁 Analyze a dataset", "🔬 Analyze a single molecule"])


# ========== SINGLE MOLECULE MODE ==========
# test smiles:
# O=C([O-])CCn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(C=NOCc4ccccc4)c3)cc21
# ========== SINGLE MOLECULE MODE ==========
if mode == "🔬 Analyze a single molecule":
    input_method = st.radio("Choose input method:", ["SMILES", "Draw Molecule"])
    single_mol = None

    if input_method == "SMILES":
        smiles = st.text_input("Enter a SMILES string:")
        if smiles:
            single_mol = Chem.MolFromSmiles(smiles)
    else:
        mol_block = st.text_area("Paste MOL block of the molecule:")
        if mol_block:
            single_mol = Chem.MolFromMolBlock(mol_block)

    if single_mol:
        st.success("Molecule loaded successfully!")

        num_confs = st.slider("Number of conformers to generate:", min_value=5, max_value=300, value=50)

        mol_H = Chem.AddHs(single_mol)
        cids = AllChem.EmbedMultipleConfs(mol_H, numConfs=num_confs, randomSeed=42)
        _ = AllChem.UFFOptimizeMoleculeConfs(mol_H)

        # Compute RMSD distance matrix
        dists = []
        for i in range(len(cids)):
            for j in range(i):
                rms = rdMolAlign.GetBestRMS(mol_H, mol_H, i, j)
                dists.append(rms)

        # Perform Butina clustering
        from rdkit.ML.Cluster import Butina
        clusters = Butina.ClusterData(dists, len(cids), 1.5, isDistData=True, reordering=True)

        st.success(f"🧬 Found {len(clusters)} clusters of conformers")

        st.markdown("### 🧪 Cluster Centroids")

        for i, cluster in enumerate(clusters):
            st.markdown(f"**Cluster {i+1}**")

            # Pick centroid: conformer with lowest average RMSD to others
            best_conf = cluster[0]
            if len(cluster) > 1:
                avg_rmsd = []
                for c1 in cluster:
                    rmsd_sum = sum(rdMolAlign.GetBestRMS(mol_H, mol_H, c1, c2) for c2 in cluster if c1 != c2)
                    avg_rmsd.append((c1, rmsd_sum / (len(cluster) - 1)))
                best_conf = min(avg_rmsd, key=lambda x: x[1])[0]

            # Align all conformers in cluster to centroid
            for c in cluster:
                rdMolAlign.AlignMol(mol_H, mol_H, prbCid=c, refCid=best_conf)

            # Visualize centroid conformer
            viewer = py3Dmol.view(width=300, height=300)
            mb = Chem.MolToMolBlock(mol_H, confId=best_conf)
            viewer.addModel(mb, "mol")
            viewer.setStyle({'stick': {}})
            viewer.setBackgroundColor("white")
            viewer.zoomTo()
            st.components.v1.html(viewer._make_html(), height=300)

# ========== DATASET MODE ==================================================================================
else:
    uploaded_file = st.file_uploader("Upload a molecule file (.sdf, .mol, .csv with SMILES)", type=["sdf", "mol", "csv"])
    mols, smiles_list = [], []

    if uploaded_file:
        ext = uploaded_file.name.split(".")[-1].lower()
        if ext == "csv":
            df = pd.read_csv(uploaded_file)
            smiles_col = next((col for col in df.columns if col.lower() in ["smiles", "smile"]), None)
            if smiles_col:
                smiles_list = df[smiles_col].dropna().tolist()
                mols, smiles_list = clean_smiles_list(smiles_list)
            else:
                st.error("No SMILES column found!")
        elif ext == "sdf":
            suppl = Chem.ForwardSDMolSupplier(uploaded_file)
            mols = [m for m in suppl if m]
            smiles_list = [Chem.MolToSmiles(m) for m in mols]
        elif ext == "mol":
            mol = Chem.MolFromMolBlock(uploaded_file.read().decode("utf-8"))
            if mol:
                mols = [mol]
                smiles_list = [Chem.MolToSmiles(mol)]

    if mols:
        with st.spinner("🔄 Analyzing molecules..."):
            if st.toggle("🔹 Use only 1000 molecules (for speed)"):
                mols = mols[:1000]
                smiles_list = smiles_list[:1000]

            fps = [get_fingerprint(m) for m in mols]
            n = len(fps)
            sim_matrix = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            dist_matrix = 1 - sim_matrix

            coords = PCA(n_components=min(10, dist_matrix.shape[1])).fit_transform(dist_matrix)

            best_score = -1
            best_k = 2
            for k in range(2, min(11, n)):
                kmeans = KMeans(n_clusters=k, random_state=0).fit(coords)
                score = silhouette_score(coords, kmeans.labels_)
                if score > best_score:
                    best_score = score
                    best_k = k

            st.success(f"✅ Found {best_k} clusters based on Silhouette Score")

            model = KMeans(n_clusters=best_k, random_state=0).fit(coords)
            labels = model.labels_

            pca_df = pd.DataFrame(coords[:, :2], columns=["PCA1", "PCA2"])
            pca_df["Cluster"] = labels
            pca_df["SMILES"] = smiles_list

            fig = px.scatter(pca_df, x="PCA1", y="PCA2", color=pca_df["Cluster"].astype(str),
                             hover_data=["SMILES"], title="Molecule Clusters",
                             color_discrete_sequence=px.colors.qualitative.Set3)
            fig.update_layout(title={"x": 0.5, "font": {"size": 24}}, plot_bgcolor="#ffffff", paper_bgcolor="#ffffff")
            st.plotly_chart(fig, use_container_width=True)

            cluster_props_summary = {}
            for clust in sorted(set(labels)):
                indices = [i for i, x in enumerate(labels) if x == clust]
                props = [calculate_properties(mols[i], smiles_list[i]) for i in indices]
                cluster_props_summary[clust] = pd.DataFrame(props).select_dtypes(include=[np.number]).mean()

            st.markdown("### 🔍 Quick Cluster Search by Properties")
            selected_props = st.multiselect("Select properties you want high values for:",
                                            ["Molecular Weight", "LogP", "H-Bond Donors", "H-Bond Acceptors",
                                             "TPSA", "Rotatable Bonds", "Aromatic Rings"])

            if selected_props:
                filtered_clusters = []
                for c, summary in cluster_props_summary.items():
                    if all(summary[prop] >= np.nanmean([cluster_props_summary[c][prop] for c in cluster_props_summary]) for prop in selected_props):
                        filtered_clusters.append(int(c))
                if filtered_clusters:
                    st.info(f"Suggested Clusters matching your properties: {filtered_clusters}")
                else:
                    st.warning("No matching cluster found.")

            selected_cluster = st.selectbox("Or select a Cluster to Explore:", sorted(pca_df["Cluster"].unique()))

            cluster_indices = pca_df[pca_df["Cluster"] == selected_cluster].index.tolist()
            st.success(f"🔍 Found {len(cluster_indices)} molecules in Cluster {selected_cluster}")

            cluster_props = []
            for idx in cluster_indices:
                mol = mols[idx]
                smile = smiles_list[idx]
                props = calculate_properties(mol, mol_name=smile)
                cluster_props.append(props)

            cluster_df = pd.DataFrame(cluster_props)
            avg_mw = cluster_df["Molecular Weight"].mean()
            avg_logp = cluster_df["LogP"].mean()
            avg_donors = cluster_df["H-Bond Donors"].mean()
            avg_acceptors = cluster_df["H-Bond Acceptors"].mean()

            size_desc = "small" if avg_mw < 300 else "medium-sized" if avg_mw < 500 else "large"
            polarity_desc = "hydrophobic (lipophilic)" if avg_logp > 2 else "hydrophilic (polar)" if avg_logp < 0 else "moderately polar"

            st.markdown(f"""<div style='padding:10px; background-color:white; border-radius:10px;'>
            <b>This cluster contains {size_desc} molecules that are {polarity_desc}.</b><br><br>
            • <b>Average Molecular Weight:</b> {avg_mw:.1f} g/mol<br>
            • <b>Average LogP:</b> {avg_logp:.2f}<br>
            • <b>Average H-Bond Donors:</b> {avg_donors:.1f}<br>
            • <b>Average H-Bond Acceptors:</b> {avg_acceptors:.1f}
            </div>""", unsafe_allow_html=True)

            for idx in cluster_indices:
                mol = mols[idx]
                smile = smiles_list[idx]
                props = calculate_properties(mol, mol_name=smile)
                st.markdown(f"<h2>Molecule {idx+1}</h2>", unsafe_allow_html=True)
                col1, col2 = st.columns([1, 2])
                with col1:
                    st.image(Draw.MolToImage(mol, size=(300, 300)))
                    viewer = show_3d_molecule(mol)
                    st.components.v1.html(viewer._make_html(), height=250)
                with col2:
                    st.dataframe(pd.DataFrame(props.items(), columns=["Property", "Value"]))

            csv = cluster_df.to_csv(index=False).encode('utf-8')
            st.download_button("Download Cluster Molecules", data=csv,
                               file_name=f"cluster_{selected_cluster}_molecules.csv", mime="text/csv")
