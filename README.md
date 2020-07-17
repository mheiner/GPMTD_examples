# GPMTD Example Code

This folder contains code to fit the Gaussian process MTD model developed in the associated article by Heiner and Kottas. All MCMC is run in [Julia](https://julialang.org/) v1.1+. One can work directly with the .jl scripts, or call them with the .sh scripts.

## Packages

Before running the .jl scripts, it is necessary to install the packages that appear after all instances of `using` in the .jl scripts. To do this, open the Julia REPL and enter package mode by pressing `]`. In package mode, enter `add` followed by each package name, separated by spaces. The *BayesInference*, *SparseProbVec*, and *MTD* packages can be installed with

```julia
pkg> add "https://github.com/mheiner/BayesInference.jl"

pkg> add "https://github.com/mheiner/SparseProbVec.jl"

pkg> add "https://github.com/mheiner/MTD.jl"
```
Other required packages include: *BSON*, *Dates*, *Distributed*, *PDMats*, *Plotly*, *Printf*, *Random*, *Random123*, *RCall*, *Statistics*, and *StatsBase*.

Post processing is run with R scripts, which require the *coda* package.

## Data

The data folder contains the two simulations described in the paper as well as the pink salmon data set (in the public domain, see <https://inport.nmfs.noaa.gov/inport/item/17256>).

## Model

The .jl scripts contain options for prior and MCMC settings/initialization. Models can be run interactively by un-commenting appropriate lines in the first block and ignoring the `ARGS` assignments.

Much of post-processing is computed with the `postProcess_loop.R` script, which can be run by sourcing the `runpostproc.sh` shell script with a trailing command-line argument identifying the posterior simulations desired with a regular expression. Results are saved in the plots folder.

Interactive post-processing can be run with `postProcess.jl`.

## Folders

- data: Contains data and scripts for simulation.
- logs: Shell output from running .jl scripts can be collected in this folder.
- plots: R post processing saves plots in this folder.
- postsim: Posterior simulations are saved to this folder.
- postsimB: Posterior simulations are moved to this folder by `postsim_to_R.jl`.
- postsim_progress: Files reporting MCMC progress and statistics are saved in this folder.
