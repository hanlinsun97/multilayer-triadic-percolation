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

@misc{sun2025triadicpercolationmultilayernetworks,
      title={Triadic percolation on multilayer networks}, 
      author={Hanlin Sun and Filippo Radicchi and Ginestra Bianconi},
      year={2025},
      eprint={2510.09341},
      archivePrefix={arXiv},
      primaryClass={nlin.AO},
      url={https://arxiv.org/abs/2510.09341}, 
}
```

