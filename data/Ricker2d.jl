
function trans!(state, r, a, b)
    x = state[1]
    y = state[2]
    x_out = x * exp( r - a*x - b*y )
    y_out = y * exp( r - a*y + b*x )

    [x_out, y_out]
end

function simulate(inits, r, a, b, nburn=Int(1e3), nkeep=Int(50e3))

    state = deepcopy(inits)

    for i = 1:nburn
        state = trans!(state, r, a, b)
    end
    println(state)

    sim = zeros(Float64, nkeep, 2)

    for i = 1:nkeep
        state = trans!(state, r, a, b)
        sim[i,:] = copy(state)
    end

    sim
end

using Random
Random.seed!(1)

## dynamics are deterministic
sim = simulate([5.0, 5.0], 2.75, 0.5, 0.07, 1000, 10010)

using Plotly

nplot = 1000

plot(sim[1:nplot,1])
plot(sim[1:nplot,1], sim[2:(nplot+1),1], mode="markers")

trace = scatter3d(Dict(
    [ :x => sim[1:(nplot),1], :y => sim[1:(nplot),2], :z => sim[2:(nplot+1),1],
      :mode => "markers",
      :marker => Dict(
        :size => 2
      ) ]
))

Plotly.plot( trace )


trace2 = scatter3d(Dict(
    [ :x => sim[1:(nplot),1], :y => sim[2:(nplot+1),1], :z => sim[3:(nplot+2),1],
      :mode => "markers",
      :marker => Dict(
        :size => 2
      ) ]
))

Plotly.plot( trace2 )


using BayesInference
x = embed(sim[:,1], 10)

y = deepcopy(x[:,1])
X = deepcopy(x[:,2:11])

using RCall
@rput y X
R"X = log(X); y = log(y); save(file='Ricker2d_log.rda', X, y)"

using BSON
bson("Ricker2d.bson", y=deepcopy(y), X=deepcopy(X))
bson("Ricker2d_log.bson", y=deepcopy(log.(y)), X=deepcopy(log.(X)))


### plots
using Plotly

r = 2.75
a = 0.5
b = 0.07

nplot = 3000

trace_sysdat = scatter3d(Dict(
    :x => (sim[1:nplot,1]),
    :y => (sim[1:nplot,2]),
    :z => (sim[2:(nplot+1),1]),
    :marker => Dict(
        :color => "rgba(0,0,0,0.70)",
        :size => 2
    ),
    :mode => "markers"
))

trace_logsysdat = scatter3d(Dict(
    :x => log.(sim[1:nplot,1]),
    :y => log.(sim[1:nplot,2]),
    :z => log.(sim[2:(nplot+1),1]),
    :marker => Dict(
        :color => "rgba(0,0,0,0.70)",
        :size => 2
    ),
    :mode => "markers"
))

trace_tdedat = scatter3d(Dict(
    :x => (X[1:nplot,1]),
    :y => (X[1:nplot,2]),
    :z => (y[1:nplot]),
    :marker => Dict(
        :color => "rgba(0,0,0,0.70)",
        :size => 2
    ),
    :mode => "markers"
))

trace_logtdedat = scatter3d(Dict(
    :x => log.(X[1:nplot,1]),
    :y => log.(X[1:nplot,2]),
    :z => log.(y[1:nplot]),
    :marker => Dict(
        :color => "rgba(0,0,0,0.70)",
        :size => 2
    ),
    :mode => "markers"
))


function sys_surf(x, y, r, a, b)
    x_out = x * exp( r - a*x - b*y )
    y_out = y * exp( r - a*y + b*x )

    [x_out, y_out]
end

function sys_surflog(lx, ly, r, a, b)
    x = exp(lx)
    y = exp(ly)
    lx_out = lx + r - a*x - b*y
    ly_out = ly + r - a*y + b*x

    [lx_out, ly_out]
end

function tde_surf(x1, x2, r, a, b)
    lx1 = log(x1)
    lx2 = log(x2)
    ax2rl = a*x2 - r + lx1 - lx2
    x1 * exp( r - a*x1 + ax2rl * exp( r + a*ax2rl/b + b*x2 ) )
end

function tde_surflog(lx1, lx2, r, a, b)
    x1 = exp(lx1)
    x2 = exp(lx2)
    ax2rl = a*x2 - r + lx1 - lx2
    lx1 + r - a*x1 + ax2rl * exp( r + a*ax2rl/b + b*x2 )
end

n_star = 900
x_star = collect(range(0.001, length=n_star, stop=(13.0)))
y_star = collect(range(0.001, length=n_star, stop=(30.0)))
lx_star = collect(range(-2.0, length=n_star, stop=(3.25)))
ly_star = collect(range(-8.0, length=n_star, stop=(3.5)))

sys_mat = zeros(Float64, n_star, n_star)
logsys_mat = zeros(Float64, n_star, n_star)
tde_mat = zeros(Float64, n_star, n_star)
logtde_mat = zeros(Float64, n_star, n_star)
for i = 1:n_star
    for j = 1:n_star
        xx, tmp = sys_surf(x_star[i], y_star[j], r, a, b)
        global sys_mat[i,j] = isinf(xx) || abs(xx) > 15.0 ? NaN : xx

        xx, tmp = sys_surflog(lx_star[i], ly_star[j], r, a, b)
        global logsys_mat[i,j] = isinf(xx) || abs(xx) > 10.0 ? NaN : xx

        xx = tde_surf(x_star[i], x_star[j], r, a, b)
        global tde_mat[i,j] = isinf(xx) || abs(xx) > 15.0 ? NaN : xx

        xx = tde_surflog(lx_star[i], lx_star[j], r, a, b)
        global logtde_mat[i,j] = isinf(xx) || abs(xx) > 3.0 ? NaN : xx
    end
end


trace_sys = surface(Dict(
  :x => x_star,
  :y => y_star,
  :z => sys_mat,
  # :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))

trace_logsys = surface(Dict(
  :x => lx_star,
  :y => ly_star,
  :z => logsys_mat,
  # :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))

trace_tde = surface(Dict(
  :x => x_star,
  :y => x_star,
  :z => tde_mat,
  # :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))

trace_logtde = surface(Dict(
  :x => lx_star,
  :y => lx_star,
  :z => logtde_mat,
  # :colorscale => "Viridis",
  :opacity => 0.7,
  :showscale => false,
  :type => "surface"
))


layout = Layout(Dict(
    :hovermode => "closest",  :margin => Dict( :r => 5, :t => 5, :b => 5, :l => 5),
    :scene => Dict(
        :xaxis => Dict(:title => "y[t-1]"),
        :yaxis => Dict(:title => "z[t-1]"),
        :zaxis => Dict(:title => "y[t]") )
    ))
data = [trace_sysdat, trace_sys]
plot(data, layout)


layout = Layout(Dict(
    :hovermode => "closest",
    :scene => Dict(
        :xaxis => Dict(:title => "log y[t-1]"),
        :yaxis => Dict(:title => "log z[t-1]"),
        :zaxis => Dict(:title => "log y[t]") )
    ))
data = [trace_logsysdat, trace_logsys]
plot(data, layout)


layout = Layout(Dict(
    :hovermode => "closest",
    :scene => Dict(
        :xaxis => Dict(:title => "y[t-1]"),
        :yaxis => Dict(:title => "y[t-2]"),
        :zaxis => Dict(:title => "y[t]") )
    ))
data = [trace_tdedat, trace_tde]
plot(data, layout)


layout = Layout(Dict(
    :hovermode => "closest",
    :scene => Dict(
        :xaxis => Dict(:title => "log y[t-1]"),
        :yaxis => Dict(:title => "log y[t-2]"),
        :zaxis => Dict(:title => "log y[t]") )
    ))
data = [trace_logtdedat, trace_logtde]
plot(data, layout)
