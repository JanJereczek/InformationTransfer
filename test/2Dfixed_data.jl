include("../scripts/intro.jl")
using DifferentialEquations
using CairoMakie
using DelimitedFiles

n_bootstrap = 5
X = readdlm("test_timeseries.csv", ',', Float64, '\n')
F = readdlm("test_F.csv", ',', Float64, '\n')
nvar, nt = size(X)
dt_out = 0.1

# Compute time derivative of states
dXdt = forward_euler(X, 1, dt_out)
T, tau, R, error_T, error_tau, error_R = bootstrapped_lianginfo_transfer(
    X, dt_out, n_bootstrap)
display(T)
display(tau)
display(F)