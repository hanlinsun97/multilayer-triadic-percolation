# Triadic percolation on multilayer networks

Julia codes for **multilayer triadic percolation** with interlayer/intralayer and positive/negative regulations written by Hanlin Sun (<hanlinsun.work@gmail.com>).

This repo contains:
- `multilayer_triadic_MCMC.jl` — Monte Carlo simulation on generated multilayer networks.
- `multilayer_triadic_theory.jl` — Mean-field theory for the corresponding coupled self-consistent dynamics.

---

We consider two structural layers **A** and **B** with:
- Poisson structural graphs with average degrees `c_A`, `c_B`
- regulatory layers (positive/negative) that can be:
  - **intra-layer**: with averaged regulatory degree `c_AA` and `c_BB`.
  - **inter-layer**: with averaged regulatory degree `c_AB` and `c_BA`.

The output is a time series of the relative size of the giant components in layer A and B `(R_A(t), R_B(t))` and orbit diagram.



---

## Requirements

**Julia packages used**
- `Plots`
- `Random`
- `DelimitedFiles`
- `LinearAlgebra` (theory code)
- `LaTeXStrings`
- `Base.Threads` (Monte Carlo code uses `@threads`)

---

# Citing
If you find the codes useful in your research, please cite the following paper:

```latex

@article{Sun2023,
author={Sun, Hanlin and Radicchi, Filippo and Kurths, J{\"u}rgen and Bianconi, Ginestra},
title={The dynamic nature of percolation on networks with triadic interactions},
journal={Nature Communications},
year={2023},
month={Mar},
day={10},
volume={14},
number={1},
pages={1308},
issn={2041-1723},
doi={10.1038/s41467-023-37019-5},
url={https://doi.org/10.1038/s41467-023-37019-5}
}
```

# License
This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).
See the LICENSE file for details.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
