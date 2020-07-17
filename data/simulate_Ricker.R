rm(list=ls())

f = function(x, theta) x * exp( theta - x)
curve(f(x, 2.6), from=-0.2, to=7)

L = 10

TT = 10e3
nburn = 100
ninit = 5
n = ninit + nburn + TT + L

sig = 0.09
theta = 2.6




set.seed(3)

y0 = numeric(n)
y0[1:ninit] = exp(rnorm(ninit, 0.0, sd=sig))

for (t in (ninit+1):n) {
	y0[t] = f(y0[t-2], theta) + rnorm(1, 0.0, sd=sig)
}

y1 = y0[-c(1:(ninit+nburn))]

length(y1) == (TT + L)

yX = embed(y1, L+1)

y = yX[,1]
X = yX[,-1]

str(y)
str(X)

# nplot = 200
# plot.ts(y[1:nplot])

# par(mfrow=c(3,2))
# for (i in 1:5) plot(X[1:nplot,i], y[1:nplot], main=paste("Lag", i))

# library("rgl")
# plot3d(X[1:nplot,3], X[1:nplot,2], X[1:nplot,1], type='p', xlab="lag 2", ylab="lag 1", zlab="y")
# dev.off()

## for validation

nprime = 1000
tprime = sort(sample(1001:10e3, nprime, replace=FALSE))

y_valid = y[tprime]
X_valid = X[tprime,]

nsim_valid = 2000

Y_valid = matrix(NA, nrow=nsim_valid, ncol=nprime)
for (i in 1:nprime) {
	Y_valid[,i] = f(X_valid[i,2], theta) + rnorm(nsim_valid, 0.0, sd=sig)
}


save(file="Ricker_single.rda", y, X, y_valid, X_valid, Y_valid, tprime)



rm(list=ls())
