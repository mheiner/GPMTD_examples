using GPMTD
using Distributed

# using BSON: @save, @load
using BSON

using Random
using PDMats
using LinearAlgebra
BLAS.set_num_threads(1)

using Statistics
using Dates

const srnd = parse(Int64, ARGS[1])
const n = parse(Int64, ARGS[2])
const L = parse(Int64, ARGS[3])
const whichdat = ARGS[4]

## Keyword args
# srnd = 427191
# n = 50
# L = 2

## Settings in script
n_burn = Int(50e3)
n_keep = 1000
n_adapt = 2000
n_thin = 20
centerscale = false

## Load data
# whichdat = "Ricker_lag2only"

BSON.@load "data/$(whichdat).bson" y X

y = deepcopy(y[1:n])
X = deepcopy(X[1:n,1:L])

## Setup
mesg = "$(whichdat)_n$(n)_L$(L)_$(srnd)"
Random.seed!(srnd)

mean(y)
std(y)

if centerscale
    X = (X .- mean(y)) ./ std(y)
    y = (y .- mean(y)) ./ std(y)
end

## Prior
prior = Prior_GPMTD(L, InterceptNormal(), MixComponentNormal())
# for j = 1:L # "deterministic"
#     prior.mixcomps[j].σ2.ν = float(n)
#     prior.mixcomps[j].σ2.s0 = 0.0001
# end

## Init
init = State_GPMTD(X, InterceptNormal(), MixComponentNormal())

## log file
logfilename = "postsim_progress/out_prog_$(mesg).txt"

## initialize model
model = Model_GPMTD(y, X, prior, init)

monitor = [:intercept, :lλ, :ζ, :μ, :σ2, :κ, :corParams, :κ_hypers, :corHypers, :fx]
# monitor = [:intercept, :lλ, :ζ :μ, :σ2, :κ, :corParams, :κ_hypers, :corHypers, :fx]

updatevars = [:intercept, :lλ, :ζ, :μ, :σ2, :κ, :corParams, :κ_hypers, :corHypers, :fx]
# updatevars = [:intercept, :lλ, :ζ, :μ, :σ2, :κ, :corParams, :κ_hypers, :corHypers, :fx]

## set up cores for distributed updates of mixcomps
# n_workers = min(length(Sys.cpu_info()) - 1, L)
n_workers = 0
addprocs(n_workers)
n_procs = nprocs()
# n_procs = 1
workers()
@everywhere using Random123, GPMTD

## runs
@time iter, accptr = mcmc!(model, 1000, save=false, n_procs=n_procs,
        report_filename=logfilename, report_freq=100, update=updatevars)
println(model.state.iter)
println(accptr)

adapt!(model, n_iter_collectSS=n_adapt, n_iter_scale=n_adapt,
    adjust_bnds=[0.01, 5.0], maxtries=10;
    report_filename=logfilename, update=updatevars)

## burn-in
timestartburn = Dates.now()
iter, accptr = mcmc!(model, n_burn, save=false, n_procs=n_procs,
    report_filename=logfilename, report_freq=10000, update=updatevars)
etr(timestartburn, n_keep, n_thin, logfilename)

sims, accptr = mcmc!(model, n_keep, save=true, thin=n_thin,
    n_procs=n_procs,
    report_filename=logfilename, report_freq=10000,
    update=updatevars, monitor=monitor)

if n_workers > 0
    rmprocs(workers())
end

using BSON
bson("postsim/postsim_$(mesg).bson",
    model=deepcopy(model), sims=deepcopy(sims),
    init=deepcopy(init), mesg=deepcopy(mesg) )
