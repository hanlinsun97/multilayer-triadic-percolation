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

@article{yvtg-wnn4,
  title = {Triadic percolation on multilayer networks},
  author = {Sun, Hanlin and Radicchi, Filippo and Bianconi, Ginestra},
  journal = {Phys. Rev. E},
  volume = {113},
  issue = {1},
  pages = {014313},
  numpages = {13},
  year = {2026},
  month = {Jan},
  publisher = {American Physical Society},
  doi = {10.1103/yvtg-wnn4},
  url = {https://link.aps.org/doi/10.1103/yvtg-wnn4}
}

```

