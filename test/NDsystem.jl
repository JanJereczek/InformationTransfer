include("../scripts/intro.jl")
using DifferentialEquations
using CairoMakie
using DelimitedFiles

function coupling_matrix(nvar::Int)
    sigma = diagm( -rand(nvar) )    # eigenvalues ∈ [-1, 0] for stable system
    Q = rand(nvar, nvar)            # random eigenvectors
    return Q * sigma * inv(Q)
end

function noise_vector(nvar::Int)
    return rand(nvar)
end

function f!(du, u, p, t)
    du[:] .= p.F * u
end

function g!(du, u, p, t)
    du[:] .= p.g
end

nvar = 2
n_bootstrap = 5
F = coupling_matrix(nvar)
g = noise_vector(nvar)
p = (F = F, g = g)

tspan = (0.0, 10.0)
dt_out = 0.1
load_data = false

u0 = zeros(nvar)
prob = SDEProblem(f!, g!, u0, tspan, p)
sol = solve(prob, saveat=dt_out)
X = hcat(sol.u...)
writedlm("test_timeseries.csv", X, ',')
writedlm("test_F.csv", F, ',')

fig, ax, l1 = lines(X[1, :])
l2 = lines!(X[2, :])
fig

# Compute time derivative of states
dXdt = forward_euler(X, 1, dt_out)
# Compute transfer information for original time series
T, tau, R = lianginfo_transfer(X, dXdt, dt)
display(T)
display(F)

T, tau, R, error_T, error_tau, error_R = bootstrapped_lianginfo_transfer(
    X, dt_out, n_bootstrap)