include("../scripts/intro.jl")
using DelimitedFiles

X = readdlm("test_timeseries.csv", ',', Float64, '\n')
nvar, nt = size(X)
dXdt = forward_euler(X, 1, dt_out)


dC = zeros(nvar, nvar)
for i in 1:nvar, j in 1:nvar
    dC[j, i] = sum( (X[j, :] .- mean(X[j, :])) .*
        (dXdt[i, :] .- mean(dXdt[i, :])) ) / (nt-1)
end

dC ≈ cov(X', dXdt')