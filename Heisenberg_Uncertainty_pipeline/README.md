SPDX-License-Identifier: MIT

Copyright (c) 2025 Stefan Len

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17356922.svg)](https://doi.org/10.5281/zenodo.17356922)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![arXiv](https://img.shields.io/badge/arXiv-physics.comp--ph-b31b1b.svg)](https://arxiv.org/)

# Heisenberg Uncertainty — Minimal Simulation

**Author:** Stefan Len  

---

## Overview
This project provides a **minimal numerical demonstration** of the **Heisenberg Uncertainty Principle**.  
It uses Gaussian wave packets to compute position- and momentum-space spreads (Δx, Δp) and verifies that:

$$
\Delta x \cdot \Delta p \; \geq \; \hbar/2
$$

The code also includes an optional **free-particle evolution** to show wave packet spreading.

---

## Features
- Generate Gaussian wave packets in position space  
- Compute uncertainties Δx and Δp via FFT  
- Sweep over initial spreads σₓ and output CSV data  
- Verify the uncertainty relation graphically  
- Optional time evolution using a split-operator method  

---

## Requirements
- Python 3.8+  
- `numpy`, `matplotlib`  

Install dependencies:
```bash
pip install numpy matplotlib
```

## License

This project is licensed under the MIT License – see the [LICENSE](../LICENSE) file for details.

## Contact

Got questions, ideas, or feedback?  
Drop me an email at **tqe.simulation@gmail.com** 
