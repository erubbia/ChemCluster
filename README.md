<p align="center">
  <img width="450" alt="Logo ChemCluster" src="https://raw.githubusercontent.com/Romainguich/ChemCluster/main/assets/Logo%20ChemCluster.png">
</p>

# - ChemCluster -

**ChemCluster** is an interactive web application for cheminformatics and molecular analysis, focusing on forming and visualizing molecular clusters built using **Streamlit**, **RDKit**, and **scikit-learn**.

Final project for the course **Practical Programming in Chemistry** — EPFL CH-200

## 📦 Package overview

**ChemCluster** is an interactive cheminformatics platform developed at **EPFL** in 2025 as part of the *Practical Programming in Chemistry* course. It is a user-friendly web application designed to explore and analyze chemical structures, either individually via the formation of conformers or as datasets. 

This tool enables users to compute key molecular properties, visualize 2D and 3D structures, and perform clustering based on molecular similarity or conformer geometry. It also offers filtering options to help select clusters matching specific physicochemical criteria.


## 🌟 Features

-  Load molecule files (`.sdf`, `.mol`) or SMILES in `.csv` format
-  Compute key molecular descriptors (MW, logP, TPSA, H-bond donors/acceptors, etc.)
-  Visualize 2D molecular structures with RDKit
-  Generate 3D conformers and visualize them interactively using Py3Dmol
-  Apply PCA for dimensionality reduction
-  Cluster molecules using KMeans with automatic silhouette score optimization
-  Click to view molecular properties directly from PCA plot
-  Export clusters and molecular data to `.csv`

## 🛠️ Installation

1. Clone the repository:

```
git clone https://github.com/erubbia/ChemCluster.git
cd ChemCluster
```

2. Create and activate the conda environment:

```
conda env create -f environment.yml
conda activate chemcluster-env
```

3. Run the Streamlit application:

```
streamlit run app.py
```

## 📖 Usage

After launching the app, access it via Streamlit’s local interface.

You can:
- Analyze a single molecule by inputting a SMILES string or drawing the structure
- Upload a dataset of molecules to perform PCA and clustering
- Click on any point in the scatter plot to view its structure and properties
- Use filters to identify clusters with desirable properties (e.g., high LogP, low MW)
- Export selected clusters as CSV files for further analysis

## 📂 License

[MIT License](LICENSE)


---

### 👨‍🔬 Developers

- Elisa Rubbia, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - erubbia](https://img.shields.io/badge/GitHub-erubbia-181717.svg?style=flat&logo=github)](https://github.com/erubbia)

- Romain Guichonnet, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - Romainguich](https://img.shields.io/badge/GitHub-Romainguich-181717.svg?style=flat&logo=github)](https://github.com/Romainguich)

- Flavia Zabala Perez, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - Flaviazab](https://img.shields.io/badge/GitHub-Flaviazab-181717.svg?style=flat&logo=github)](https://github.com/Flaviazab)
