using BSON
using RCall

R"load('Ricker_single.rda')"
R"ls()"

@rget y X y_valid X_valid Y_valid tprime

bson("Ricker_single.bson", y = deepcopy(y),
    X = deepcopy(X), y_valid = deepcopy(y_valid),
    X_valid = deepcopy(X_valid), Y_valid = deepcopy(Y_valid),
    tprime = deepcopy(tprime))
