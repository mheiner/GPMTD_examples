using GPMTD
using LinearAlgebra
using SparseProbVec
using PDMats
using Random123
using Statistics

using BSON
using BSON: @load
whichdat = "postsim_Ricker_lag2only_n50_L2_427191"
BSON.@load "postsim/$(whichdat).bson" model sims init mesg
# BSON.@load "postsimB/$(whichdat).bson" model sims init mesg

using Plotly
using StatsBase

nsim = length(sims)

iter_use = 1:nsim
iter_use = range(1, nsim, step=2)
iter_use = 1900:nsim

# counts of mixture component membership
sims_Scounts = permutedims( hcat( [ StatsBase.counts( sims[i][:ζ], 0:model.L ) for i in iter_use ]... ) )
mean(sims_Scounts ./ model.n, dims=1)

tmp = hcat( [ sims[ii][:ζ] for ii in iter_use ]... )
sims_Scount_obs = permutedims(hcat([ StatsBase.counts( tmp[i,:], 0:model.L ) for i = 1:model.n ]...) ./ length(iter_use))

# log-likelihood
sims_llik = [ sims[i][:llik] for i in iter_use ]
plot(sims_llik)

## parameter inferences
j = 1 # 0 for intercept
lam = [ exp(sims[i][:lλ][j+1]) for i in iter_use ]
plot(lam)

for j = 0:model.L
    lam = [ exp(sims[i][:lλ][j+1]) for i in iter_use ]
    println("lam $(j) mean ", round(mean(lam), digits=3))
    println("lam $(j) q025 ", round(quantile(lam, 0.025), digits=4))
    println("lam $(j) q975 ", round(quantile(lam, 0.975), digits=4), "\n")
end

int_mu = [ sims[i][:intercept].μ for i in iter_use ]
plot(int_mu)

int_sig2 = [ sims[i][:intercept].σ2 for i in iter_use ]
plot(int_sig2)

j = 1
mixc_mu = [ sims[i][:mixcomps][j][:μ] for i in iter_use ]
mean(mixc_mu)
std(mixc_mu)
plot(mixc_mu)

mixc_sig2 = [ sims[i][:mixcomps][j][:σ2] for i in iter_use ]
plot(mixc_sig2)
mean(mixc_sig2)
plot(sqrt.(mixc_sig2))
mean(sqrt.(mixc_sig2))

j = 1
mixc_kap = [ sims[i][:mixcomps][j][:κ] for i in iter_use ]
plot(mixc_kap)

mixc_ls = [ sims[i][:mixcomps][j][:corParams].lenscale for i in iter_use ]
plot(mixc_ls)

covstat = [ sims[i][:mixcomps][j][:κ]*sims[i][:mixcomps][j][:σ2]/(sims[i][:mixcomps][j][:corParams].lenscale^2.5) for i in iter_use ]
plot(covstat)

j = 2
kap_df = [ sims[i][:mixcomps][j][:κ_hypers].κ_ν for i in iter_use ]
plot(kap_df)

j = 2
kap0 = [ sims[i][:mixcomps][j][:κ_hypers].κ0 for i in iter_use ]
plot(kap0)

j = 2
ls_df = [ sims[i][:mixcomps][j][:corHypers].lenscale_ν for i in iter_use ]
plot(ls_df)

j = 2
ls0 = [ sims[i][:mixcomps][j][:corHypers].lenscale0 for i in iter_use ]
plot(ls0)


j = 2
ord = sortperm(model.X[:,j])
F = hcat([ sims[i][:mixcomps][j][:fx] .+ sims[i][:mixcomps][j][:μ] for i in iter_use ]...)
Fmean = mean(F, dims=2)[:,1]
Fq025 = [ quantile(F[i,:], 0.025) for i = 1:size(F,1)]
Fq975 = [ quantile(F[i,:], 0.975) for i = 1:size(F,1)]

trace_mean = scatter(Dict("x" => model.X[ord,j], "y" => Fmean[ord], "name" => "mean",
    "mode" => "both", "type" => "scatter", "marker" => Dict("color" => "black", "line" => Dict("color" => "black"))))
trace_05 = scatter(Dict("x" => model.X[ord,j], "y" => Fq025[ord], "name" => "Q025",
    "mode" => "both", "type" => "scatter", "marker" => Dict("color" => "red", "line" => Dict("color" => "red"))))
trace_95 = scatter(Dict("x" => model.X[ord,j], "y" => Fq975[ord], "name" => "Q975",
    "mode" => "both", "type" => "scatter", "marker" => Dict("color" => "orange", "line" => Dict("color" => "orange"))))

plot([trace_mean, trace_05, trace_95])

plot(scatter(Dict("x" => X[ord,j], "y" => F[ord,50],
    "mode" => "both", "type" => "scatter")))




### inferences on a grid of x values

include("Rfuns.jl")

n_star = 100
minx = minimum(model.X)
maxx = maximum(model.X)
ranx = maxx - minx

X_star_vals = collect(range(minx - 0.1*ranx, length=n_star,
                        stop=(maxx + 0.1*ranx)))
# X_star = rcopy(xpndgrid([X_star_vals for k = 1:model.K]...))
X_star_vals = collect(range(2.0, length=n_star,
                        stop=11.5))


γindx = [2]

# iter_use = sort(findall([ sims[ii][:γ] == γ_use for ii in 1:length(sims)]))
iter_use = 1:nsim
# iter_use = 1900:2000
iter_use = range(1, nsim, step=2)

Kuse = length(γindx)
X_star00 = rcopy(xpndgrid([X_star_vals for k = 1:Kuse]...))

# X_star0 = zeros(Float64, size(X_star00,1), model.K)
X_star0 = fill(mean(model.y), size(X_star00,1), model.L)
# X_star0 = rand(size(X_star00,1), model.L) .* ranx .+ minx # draw uniformly over the range

X_star0[:,γindx] = deepcopy(X_star00)
X_star = deepcopy(X_star0)

# transition mean function
@time Ey = getEyPy(sims[iter_use], X_star, model.X, model.D)

mean_Ey = mean(Ey, dims=2)[:,1]

q025_Ey = [quantile(Ey[jj,:], 0.025) for jj = 1:size(Ey,1)]
q975_Ey = [quantile(Ey[jj,:], 0.975) for jj = 1:size(Ey,1)]




## one dimension
size(Ey)
plotR(X_star_vals, mean_Ey, ylab="", xlab="",
    ylim=[minimum(q025_Ey), maximum(q975_Ey)], lwd=2, typ="l")
plotR(X_star_vals, mean_Ey, ylab="", xlab="",
    ylim=[0.0, 12.3], lwd=2, typ="l")
plotR(X_star_vals, mean_Ey, ylab="", xlab="",
    ylim=[minimum(model.y)-0.05*ranx, maximum(model.y)+0.05*ranx], lwd=2, typ="l")
polygonR(vcat(X_star_vals, reverse(X_star_vals)), vcat(q025_Ey, reverse(q975_Ey)), col="gray85", border="NA")
pointsR(model.X[:,γindx[1]], model.y)

linesR(X_star_vals, mean_Ey, lwd=2)

R"title(ylab='y')"
R"title(xlab='y lag')"

R"dev.off()"


### two dimensions
# if extreme values mess up graph
# mean_Ey[findall(mean_Ey .< -1.0)] *= 0.0

# hcat(X_star[:,1], X_star[:,2], mean_Ey)
mean_Ey_mat = Matrix(reshape(mean_Ey, fill(n_star,Kuse)...))

q025_Ey_mat = Matrix(reshape(q025_Ey, fill(n_star,Kuse)...))
q975_Ey_mat = Matrix(reshape(q975_Ey, fill(n_star,Kuse)...))


trace1 = surface(Dict(
  :x => X_star_vals,
  :y => X_star_vals,
  :z => q025_Ey_mat,
  :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))

trace2 = surface(Dict(
  :x => X_star_vals,
  :y => X_star_vals,
  :z => mean_Ey_mat,
  # :colorbar => Dict(
  #   :ticklen => 2,
  #   :title => ""
  # ),
  :colorscale => "Viridis",
  :showscale => false,
  :type => "surface"
))

trace3 = surface(Dict(
  :x => X_star_vals,
  :y => X_star_vals,
  :z => q975_Ey_mat,
  :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))

trace0 = scatter3d(Dict(
    :x => model.X[:,γindx[1]],
    :y => model.X[:,γindx[2]],
    :z => model.y,
    :marker => Dict(
        :color => "rgba(0,0,0,0.70)",
        :size => 2
    ),
    :mode => "markers"
))

tracetrue = surface(Dict(
  :x => X_star_vals,
  :y => X_star_vals,
  :z => truemean_mat,
  # :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))

data = [trace2, trace0]
data = [trace1, trace2, trace3, trace0#=, tracetest=#]
data = [tracetrue, trace2, trace0]
data = [trace0]

layout = Layout(Dict(
  :hovermode => "closest",
  :margin => Dict(
    :r => 10,
    :t => 25,
    :b => 40,
    :l => 60
  ),
  :scene => Dict(
    :xaxis => Dict(:title => "log y[t-$(γindx[1])]",
                   :range => [-1.5, 3.0] #=,
                   :type => "log" #=,
                   :autorange => true =# =#),
    :yaxis => Dict(:title => "log y[t-$(γindx[2])]",
                    :range => [-1.5, 3.0] #=,
                    :type => "log" #=,
                    :autorange => true =# =#),
    :zaxis => Dict(:title => "log y[t]",
                    :range => [-1.5, 3.0] #=,
                    :type => "log" #=,
                    :autorange => true =# =#),
  ),
  :showlegend => false
))

p = plot(data, layout)

plot(data)






## transition density
n_y = 150
y_grid = collect(range(minx - 0.3*ranx, length=n_y,
                        stop=(maxx + 0.3*ranx)))

# number of distinct x locations at which to evaluate
n_star_Py = 3

# x values
X_star_Py = deepcopy(X_star[1:n_star_Py,:])
X_star_Py[1, γindx[1]] = 2.0
X_star_Py[2, γindx[1]] = 5.0
X_star_Py[3, γindx[1]] = 7.0

X_star_Py = deepcopy(X_star[1:n_star_Py,:])
X_star_Py[1, γindx[1]] = 50.0
X_star_Py[2, γindx[1]] = 66.0
X_star_Py[3, γindx[1]] = 80.0

X_star_Py = deepcopy(X_star[1:n_star_Py,:])
X_star_Py[1, γindx[1]:γindx[2]] = [10.0, 2.0]
X_star_Py[2, γindx[1]:γindx[2]] = [1.0, 1.0]
X_star_Py[3, γindx[1]:γindx[2]] = [2.0, 7.0]

X_star_Py = deepcopy(X_star[1:n_star_Py,:])
X_star_Py[1, γindx[1]:γindx[2]] = [2.3, 0.8]
X_star_Py[2, γindx[1]:γindx[2]] = [1.0, -0.5]
X_star_Py[3, γindx[1]:γindx[2]] = [1.0, 1.6]
X_star_Py[4, γindx[1]:γindx[2]] = [1.7, 0.7]
X_star_Py[5, γindx[1]:γindx[2]] = [-0.75, 2.3]


X_star_Py

# transition mean and densities
Ey0, Py = getEyPy(sims[iter_use], X_star_Py,
            model.X, model.D; densout=true, y_grid=y_grid)

mean_Py = mean(Py, dims=3)[:,:,1]
q025_Py = [quantile(Py[i,j,:], 0.025) for i = 1:n_y, j = 1:n_star_Py]
# p50_Py = [quantile(Py[i,j,:], 0.5) for i = 1:n_y, j = 1:n_star_Py]
q975_Py = [quantile(Py[i,j,:], 0.975) for i = 1:n_y, j = 1:n_star_Py]


xind = 2
X_star_Py[xind,:]

R"par(mar=c(4,4,2,1)+0.1)"
plotR(y_grid, mean_Py[:,xind],
    main="y[t-1] = $(round(X_star_Py[xind,γindx[1]], digits=2))",
    ylab="density", xlab="y[t] (minutes)", typ="l", ylim=[0.0, maximum(q975_Py[:,xind])], lty=1)
polygonR(vcat(y_grid, reverse(y_grid)), vcat(q025_Py[:,xind], reverse(q975_Py[:,xind])), col="gray90", border="NA")
linesR(y_grid, mean_Py[:,xind], lwd=2)
R"dev.off()"

R"par(mar=c(4,4,2,1)+0.1)"
plotR(y_grid, mean_Py[:,xind],
    main="y[t-1] = $(round(X_star_Py[xind,γindx[1]], digits=2))",
    ylab="density", xlab="y[t] (minutes)", typ="l", ylim=[0.0, maximum(q975_Py[:,xind])], lty=1)
polygonR(vcat(y_grid, reverse(y_grid)), vcat(q025_Py[:,xind], reverse(q975_Py[:,xind])), col="gray90", border="NA")
linesR(y_grid, mean_Py[:,xind], lwd=2)
R"dev.off()"


R"par(mar=c(4,4,2,1)+0.1)"
plotR(y_grid, mean_Py[:,xind],
    main="y[t-1] = $(round(X_star_Py[xind,γindx[1]], digits=2)), y[t-2] = $(round(X_star_Py[xind,γindx[2]], digits=2))",
    ylab="density", xlab="y[t] (minutes)", typ="l", ylim=[0.0, maximum(q975_Py[:,xind])], lty=1)
polygonR(vcat(y_grid, reverse(y_grid)), vcat(q025_Py[:,xind], reverse(q975_Py[:,xind])), col="gray90", border="NA")
linesR(y_grid, mean_Py[:,xind], lwd=2)
R"dev.off()"
