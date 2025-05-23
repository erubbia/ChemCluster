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

- Upload .sdf, .mol, or .csv files containing SMILES
- Input or draw a single molecule and generate 3D conformers
- Compute key molecular properties (MW, logP, H-bonding, etc.)
- Visualize molecules in 2D (RDKit) and interactively in 3D (Py3Dmol)
- Cluster molecules using PCA + KMeans with silhouette score optimization
- Click points on the PCA plot to inspect molecules and properties
- Overlay and compare 3D cluster centroids for conformers
- Filter clusters based on desired property profiles
- Export results and clusters as .csv files

## 🛠️ Installation

1. Install from [PyPI](https://pypi.org/project/chemcluster/):

```bash
pip install chemcluster
```

2. Run the app:

```bash
chemcluster
```

This will open the ChemCluster interface in your browser.

### To run locally from source:

```bash
git clone https://github.com/erubbia/ChemCluster.git
cd ChemCluster
conda env create -f environment.yml
conda activate chemcluster-env
(chemcluster-env) $ pip install -e .
```
**Testing** can be done with 'pytest' or 'tox':
```bash
(chemcluster-env) $ pytest
# or 
(chemcluster-env) $ tox
```

## 📖 Usage

**Single molecule mode:**
- Draw and paste SMILES to visualize and cluster conformers
- View and overlay optimized 3D centroid structures

**Data set mode:**
- Upload a SMILES data set to analyze chemical space
- Perform PCA + KMeans clustering with property-based filters
- Click to view molecules and export clusters
  

## 📂 License

[MIT License](LICENSE)


---

### 👨‍🔬 Developers

- Elisa Rubbia, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - erubbia](https://img.shields.io/badge/GitHub-erubbia-181717.svg?style=flat&logo=github)](https://github.com/erubbia)

- Romain Guichonnet, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - Romainguich](https://img.shields.io/badge/GitHub-Romainguich-181717.svg?style=flat&logo=github)](https://github.com/Romainguich)

- Flavia Zabala Perez, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - Flaviazab](https://img.shields.io/badge/GitHub-Flaviazab-181717.svg?style=flat&logo=github)](https://github.com/Flaviazab)
