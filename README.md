# Triadic percolation on multilayer networks

Julia codes for **multilayer triadic percolation** with interlayer/intralayer and positive/negative regulations written by Hanlin Sun (<hanlinsun.work@gmail.com>).

This repo contains:
- `multilayer_triadic_MCMC.jl` — Monte Carlo simulation on generated multilayer networks. :contentReference[oaicite:0]{index=0}  
- `multilayer_triadic_theory.jl` — Mean-field theory for the corresponding coupled self-consistent dynamics:contentReference[oaicite:1]{index=1}

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

