# Circuit compression for 2D quantum dynamics

## Overview

`cc2d_example` is a repo designed to illustrate the circtuit compression algorithm for the dynamics of 2d systems introduced in [[1]](https://arxiv.org/abs/2507.01883).

## Installation
When running the `TFIM_example.ipynb` example notebook for the first time, install the needed packages.
```julia
using Pkg
Pkg.activate(".")

Pkg.add(name="PauliPropagation", version="0.3.0")
Pkg.add("ReverseDiff")
Pkg.add("OptimKit")
```
