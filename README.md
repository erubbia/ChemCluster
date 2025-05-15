<img width="300" alt="Logo MolSearch" src="https://raw.githubusercontent.com/Romainguich/ChemCluster/main/assets/Logo%20ChemCluster.png">

# - ChemCluster -

**ChemCluster** is an interactive web application for cheminformatics and molecular analysis, focusing on forming and visualizing molecular clusters built using **Streamlit**, **RDKit**, and **scikit-learn**.

Final project for the course **Practical Programming in Chemistry** — EPFL CH-200

## 📦 Package overview

**ChemCluster** is an interactive cheminformatics platform developed at **EPFL** in 2025 as part of the *Practical Programming in Chemistry* course. It is a user-friendly web application designed to explore and analyze chemical structures, either individually via the formation of conformers or as datasets. 

This tool enables users to compute key molecular properties, visualize 2D and 3D structures, and perform clustering based on molecular similarity or conformer geometry. It also offers filtering options to help select clusters matching specific physicochemical criteria.

---

### 👨‍🔬 Developers

- Elisa Rubbia, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - erubbia](https://img.shields.io/badge/GitHub-erubbia-181717.svg?style=flat&logo=github)](https://github.com/erubbia)

- Romain Guichonnet, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - Romainguich](https://img.shields.io/badge/GitHub-Romainguich-181717.svg?style=flat&logo=github)](https://github.com/Romainguich)

- Flavia Zabala Perez, Master's student in Molecular and Biological Chemistry at EPFL [![GitHub - Flaviazab](https://img.shields.io/badge/GitHub-Flaviazab-181717.svg?style=flat&logo=github)](https://github.com/Flaviazab)


# copier-liac

[Copier](https://github.com/copier-org/copier) template for pure Python libraries.

_As simple as possible. No magic._

Adapted from [copier-pylib](https://github.com/astrojuanlu/copier-pylib) and [snekpack](https://github.com/cthoyt/cookiecutter-snekpack).

## Features

- [Hatch] for packaging.
- [pytest] for testing.
- [tox] for automation of test runners and other stuff.
- [Sphinx] for documentation
- [GitHub Actions] for continuous integration

If you want more automated checks, look at https://github.com/schwallergroup/copier-liac.
- [mypy] for type checks.
- [ruff] for style checks and automatic Python code formatting.
- [pre-commit] for optional automation of style checks.

## Usage

Install `copier`:

```
pip install copier
```

Generate a Python package:

```
copier copy gh:schwallergroup/copier-liac-minimal path/to/destination
```

## License

[MIT License](LICENSE)

[copier]: https://github.com/copier-org/copier/
[Hatch]: https://hatch.pypa.io/
[pytest]: https://docs.pytest.org/
[Sphinx]: http://www.sphinx-doc.org/
[tox]: https://tox.readthedocs.io/
