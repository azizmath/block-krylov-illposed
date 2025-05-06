# Block Krylov Methods for Discrete Ill-posed Problems

This repository contains MATLAB implementations supporting the paper:

**On the Block Lanczos and Block Golub–Kahan Reduction Methods Applied to Discrete Ill-posed Problems**  
*Abdulaziz Alqahtani, Silvia Gazzola, Lothar Reichel, Giuseppe Rodriguez*  
📄 *Numerical Linear Algebra with Applications*, 2021  
🔗 [DOI: 10.1002/nla.2376](https://doi.org/10.1002/nla.2376)

---

## 🔍 Overview

This work explores how block Krylov subspace methods — specifically the **Block Lanczos** and **Block Golub–Kahan** algorithms — can be used to solve large-scale discrete ill-posed problems with multiple right-hand sides. These problems arise in applications like tomography, signal processing, and inverse heat conduction.

We show how these methods perform compared to traditional techniques, and include example code so others can easily try them out and reproduce the results.

---

## 🧩 Dependency: IR Tools

This project uses test problems and solvers from the [IR Tools MATLAB package](https://github.com/jnagy1/IRtools) by Silvia Gazzola, P. C. Hansen, and James G. Nagy.

To use:
1. Download IR Tools from GitHub: https://github.com/jnagy1/IRtools
2. Add it to your MATLAB path (e.g., using `IRtools_setup`)

---

## 💻 Repository Contents

| File / Folder | Description |
|---------------|-------------|
| `Block_Lanczos_tridiagonalization.m` | Block Lanczos tridiagonalization |
| `Block_Golub_Kahan_bidiagonalization.m` | Block Golub–Kahan bidiagonalization |
| `testbgkb.m` / `testblt.m` | Scripts to test both algorithms on example problems |
| `baart.m`, `foxgood.m`, etc. | Ill-posed problem generators from IR Tools |
| `figure3blt.m`, `figure3bgkb.m` | Scripted figures to reproduce paper plots |
| `IRtools-master/` | [Optional] External dependency folder (or download separately) |

---

## ▶️ Running Examples

To run a test using Block Golub–Kahan:

```matlab
testbgkb
```

To reproduce Figure 3 of the paper:

```matlab
figure3blt
figure3bgkb
```

Make sure IR Tools is added to your MATLAB path.

---

## 📖 Citation

If you use this code, please cite:

```bibtex
@article{alqahtani2021block,
  title={On the block Lanczos and block Golub--Kahan reduction methods applied to discrete ill-posed problems},
  author={Alqahtani, Abdulaziz and Gazzola, Silvia and Reichel, Lothar and Rodriguez, Giuseppe},
  journal={Numerical Linear Algebra with Applications},
  volume={28},
  number={5},
  pages={e2376},
  year={2021},
  publisher={Wiley}
}
```

---

## 📫 Contact

For questions, please reach out:  
📧 alqahtani.math@gmail.com  
🔗 [LinkedIn](https://www.linkedin.com/in/alqahtani-math)

---

© 2021 Abdulaziz Alqahtani. Licensed under the MIT License.
