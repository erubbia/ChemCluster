
# Streamlit and more
import streamlit as st
import pandas as pd
from io import BytesIO
import base64

# RDKit tools
from rdkit import Chem
from rdkit.Chem import (
    Descriptors, Crippen, Lipinski, rdMolDescriptors,
    Draw, AllChem, DataStructs, rdDistGeom, rdMolAlign
)
from rdkit.ML.Cluster import Butina

# Analysis
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# Visualization
import matplotlib
matplotlib.use('MacOSX')  # for local plotting compatibility on macOS
import plotly.express as px
import plotly.graph_objects as go
import py3Dmol
from streamlit_plotly_events import plotly_events

from chemcluster.core import (
    calculate_properties,
    get_fingerprint,
    clean_smiles_list,
    show_3d_molecule,
    mol_to_base64_img,
)

def main():
    # App setup
    st.set_page_config(page_title="ChemCluster", layout="wide")

    st.markdown(
        """
        <div style="text-align: center; padding-bottom: 10px;">
           <img src="https://raw.githubusercontent.com/Romainguich/ChemCluster/main/assets/Logo%20ChemCluster.png" alt="ChemCluster Logo" width="400">
        </div>
        <div style="text-align: left; padding-left: 0px; padding-bottom: 10px;">
            <span style="font-size: 24px; font-weight: bold; color: #11111;">
                Welcome to <span style="color: #1f3b57;">Chem</span><span style="color: #6ea37f;">Cluster</span> !
            </span>
        </div>
        """,
        unsafe_allow_html=True
    )

    st.markdown("""
        <style>
        html, body, [class*="css"] {
            font-family: 'Nunito', sans-serif;
            font-size: 18px;
            background-color: #f7f9fb;
        }
        h1, h2, h3 {
            color: #1f4e79;
        }
        .stButton > button {
            background-color: #BCD4CC;
            color: #1f1f1f;
            border-radius: 10px;
            padding: 0.6em 1.2em;
            font-size: 18px;
            transition: all 0.2s ease-in-out;
            border: none;
        }
        .stButton > button:hover {
            background-color: #98c1b6;
            transform: scale(1.05);
            color: #000000;
        }
        div[data-baseweb="radio"] > div > label > div:first-child,
        div[data-baseweb="checkbox"] > div > div {
            background-color: #BCD4CC;
            border-color: #BCD4CC;
        }
        div[data-baseweb="radio"] > div > label > div:first-child:hover,
        div[data-baseweb="checkbox"] > div > div:hover {
            background-color: #98c1b6;
            border-color: #98c1b6;
        }
        </style>
    """, unsafe_allow_html=True)

    mode = st.radio("Choose analysis mode:", ["Analyze a dataset", "Analyze a single molecule"])

    # ========== SINGLE MOLECULE MODE ==========
    if mode == "Analyze a single molecule":
        input_method = st.radio("Choose input method:", ["SMILES", "Draw molecule"])
        single_mol = None

        if input_method == "SMILES":
            smiles = st.text_input("Enter a SMILES string:")
            if smiles:
                single_mol = Chem.MolFromSmiles(smiles)
                if single_mol is None:
                    st.error("Invalid SMILES")

        elif input_method == "Draw molecule":
            st.write("Draw your molecule below:")
            st.info(
                "⚡ Click on 'Create new document' in the editor, draw your molecule, click the 'Save to file' button and choose the SMILES format to get the corresponding string, then paste it into the 'Choose input method: SMILES' section to visualize it!"
            )
            st.components.v1.html(
                """
                <iframe src="https://partridgejiang.github.io/Kekule.js/demos/items/chemEditor/chemEditor.html"
                        width="900" height="600" style="border:none;"></iframe>
                """,
                height=650,
            )
            

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

            st.success(f"✅ Found {len(clusters)} clusters of conformers")
            st.markdown("### Cluster Centroids")

            centroid_ids = []
            for i, cluster in enumerate(clusters):
                # Pick centroid: conformer with lowest average RMSD to others
                best_conf = cluster[0]
                if len(cluster) > 1:
                    avg_rmsd = []
                    for c1 in cluster:
                        rmsd_sum = sum(rdMolAlign.GetBestRMS(mol_H, mol_H, c1, c2) for c2 in cluster if c1 != c2)
                        avg_rmsd.append((c1, rmsd_sum / (len(cluster) - 1)))
                    best_conf = min(avg_rmsd, key=lambda x: x[1])[0]
                centroid_ids.append((i + 1, best_conf))

            # Multiselect box to pick which centroids to overlay
            st.markdown("#### 🔍 Compare Centroids")
            selected = st.multiselect(
                "Select cluster centroids to overlay in 3D:",
                [f"Cluster {i}" for i, _ in centroid_ids],
                default=[f"Cluster {i}" for i, _ in centroid_ids[:1]]
            )

            viewer = py3Dmol.view(width=400, height=400)
            for i, conf_id in centroid_ids:
                if f"Cluster {i}" in selected:
                    mb = Chem.MolToMolBlock(mol_H, confId=conf_id)
                    viewer.addModel(mb, "mol")
            viewer.setStyle({'stick': {}})
            viewer.setBackgroundColor("white")
            viewer.zoomTo()
            st.components.v1.html(viewer._make_html(), height=400)
            
    # ========== DATASET MODE ==========
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

                # Create PCA dataframe
                pca_df = pd.DataFrame(coords[:, :2], columns=["PCA1", "PCA2"])
                pca_df["Cluster"] = labels
                pca_df["SMILES"] = smiles_list
                pca_df["Cluster_Label"] = pca_df["Cluster"].astype(str)

                # Fonction image mini base64
                def mol_to_base64_img(mol, size=(100, 100)):
                    img = Draw.MolToImage(mol, size=size)
                    buffer = BytesIO()
                    img.save(buffer, format="PNG")
                    buffer.seek(0)
                    return "data:image/png;base64," + base64.b64encode(buffer.read()).decode()

                pca_df["MolImg"] = [mol_to_base64_img(m) for m in mols]
                sorted_cluster_labels = sorted(pca_df["Cluster_Label"].unique(), key=lambda x: int(x))

                fig = px.scatter(
                    pca_df,
                    x="PCA1",
                    y="PCA2",
                    color="Cluster_Label",
                    hover_data={"SMILES": True, "PCA1": False, "PCA2": False},
                    custom_data=["SMILES"],
                    title="Molecule Clusters (Click to view molecule)",
                    color_discrete_sequence=px.colors.qualitative.Pastel,
                    category_orders={"Cluster_Label": sorted_cluster_labels}
                )

                fig.update_layout(title={"x": 0.5, "font": {"size": 24}}, plot_bgcolor="#ffffff", paper_bgcolor="#ffffff")

                # Interact with the plot
                selected_points = plotly_events(fig, click_event=True, hover_event=False, override_height=600)

                # Display selected molecule (no second chart)
                if selected_points:
                    idx = selected_points[0]['pointIndex']
                    mol = mols[idx]
                    smile = smiles_list[idx]
                    st.markdown("### ➡️ Molecule Selected")
                    col1, col2 = st.columns([1, 2])
                    with col1:
                        st.image(Draw.MolToImage(mol, size=(300, 300)))
                        viewer = show_3d_molecule(mol)
                        st.components.v1.html(viewer._make_html(), height=250)
                    with col2:
                        props = calculate_properties(mol, mol_name=smile)
                        st.dataframe(pd.DataFrame(props.items(), columns=["Property", "Value"]))

                # -------- Property-Based Cluster Search --------
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


if __name__ == "__main__":
    main()