
fileid = parse(Int64, ARGS[1])
simregexp = Regex("postsim.*_.*$(fileid).*.bson")

using GPMTD
using SparseProbVec

using BSON
using BSON: @load

using RCall
using StatsBase

## for compatibility with the imported types
using Random
using Random123
using PDMats
using LinearAlgebra
using Statistics
using StatsBase
using Dates

files = filter(x -> occursin(simregexp, x), readdir("./postsim"))
nfiles = length(files)

i = 1
for i in 1:nfiles
    file_now = files[i]
    println("starting file: $(file_now)")

    BSON.@load "./postsim/$(file_now)" model sims init mesg
    nsim = length(sims)

    ### Send to R
    R"rm(list=ls())"
    iter = model.state.iter
    accptr =  model.state.accpt ./ model.state.iter
    L = model.L
    n = model.n
    priorinfo = "$(model.prior)"
    initinfo = "$(init)"

    sims_llik = [ sims[i][:llik] for i=1:nsim ]
    # sims_Scounts = permutedims( hcat( [ StatsBase.counts( sims[i][:ζ], 0:L ) for i = 1:nsim ]... ) )

    sims_lam = permutedims(hcat([ exp.(sims[i][:lλ]) for i=1:nsim ]...))
    sims_int_mu = [ sims[i][:intercept].μ for i=1:nsim ]
    sims_int_sig2 = [ sims[i][:intercept].σ2 for i=1:nsim ]

    @rput mesg nsim iter accptr L n priorinfo initinfo sims_llik sims_lam sims_int_mu sims_int_sig2

    R"ls()"
    R"save.image(file=paste0('postsim/mcmc_', $(mesg), '.rda'))"

    cp("./postsim/$(file_now)", "./postsimB/$(file_now)")
    rm("./postsim/$(file_now)")

    println("completed file number: $(i) of $(nfiles) \n\n")
end
