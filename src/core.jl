using Statistics, LinearAlgebra, Random, TimeseriesSurrogates
using .Threads

"""

    InformationTransfer

Return a struct containing the results of an information transfer analysis.
"""
struct InformationTransfer
    T::Matrix
    tau::Matrix
    R::Matrix
    error_T::Matrix
    error_tau::Matrix
    error_R::Matrix
end

"""

    forward_euler(X::Matrix{T}, k::Int, dt::Real)

Forward-Euler differentiation with spacing `k` and assuming that `X` has dimension
`nvar x nt`.
"""
function forward_euler(X::Matrix, t::Vector, k::Int)
    nt, dXdt, T = init_euler_diff(X, t)
    dXdt[:, 1:nt-k] .= lagged_difference(X, k, nt) ./ lagged_difference(T, k, nt)
    return dXdt
end

"""

    backward_euler(X::Matrix{T}, k::Int, dt::Real)

Backward-Euler differentiation with spacing `k` and assuming that `X` has dimension
`nvar x nt`.
"""
function backward_euler(X::Matrix, t::AbstractVector, k::Int)
    nt, dXdt, T = init_euler_diff(X, t)
    dXdt[:, k:nt] .= lagged_difference(X, k, nt) ./ lagged_difference(T, k, nt)
    return dXdt
end

function init_euler_diff(X::Matrix, t::AbstractVector)
    nvar, nt = size(X)
    if length(t) != nt
        error("Dimensions of data matrix and time vector must coincide.")
    end
    dXdt = zeros(eltype(X), nvar, nt)
    T = ones(nvar) * t'
    return nt, dXdt, T
end

function lagged_difference(X::Matrix, k::Int, nt::Int)
    return (view(X, :, k+1:nt) - view(X, :, 1:nt-k))
end

"""

    liang_index(detC::Real, Detla_jk::Vector, C_kdi::Vector,
    C_ij::Real, C_ii::Real)

Compute the (j,i)-entry of the information-transfer matrix as defined by Liang (2021).

## Inputs
 - `detC`: the determinant of the data covariance matrix `` C ``.
 - `Detla_jk`: the `` j ``-th row of `` \\mathrm{det}(C) \\cdot C^{-1} ``.
 - `C_kdi`: the `` i ``-th column of `` \\mathrm{cov}(X, \\dfrac{\\mathrm{d}X}{\\mathrm{d}t}) ``.
 - `C_ij`: the entry ``(i,j)`` of `` C ``.
 - `C_ii`: the entry ``(i,i)`` of `` C ``.
"""
function liang_index(detC::Real, Detla_jk::Vector, C_kdi::Vector,
    C_ij::Real, C_ii::Real)
    return (1 / detC) * sum(Detla_jk .* C_kdi) * (C_ij / C_ii)
end

"""

    liang_normindex(detC::Real, Delta_ik::Vector, C_kdi::Vector,
        T_all::Vector, T_ii::Real, g_ii::Real, C_ii::Real, T_ji::Real)

Compute the (j,i)-entry of the normalized information-transfer matrix as defined by Liang (2021).
"""
function liang_normindex(detC::Real, Delta_ik::Vector, C_kdi::Vector,
    T_all::Vector, T_ii::Real, g_ii::Real, C_ii::Real, T_ji::Real)
        # self-contribution (eq. 15)
        selfcontrib = (1 / detC) * sum(Delta_ik .* C_kdi)
        # all other transfers contribution (eq. 20)
        transfer = sum(abs.(T_all)) - abs(T_ii)
        # noise contribution
        noise = 0.5 * g_ii / C_ii
        # normalizer (eq. 20)
        Z = abs(selfcontrib) + transfer + abs(noise)
        # relative rate of information flowing from xj to xi (eq. 19)
        return 100 * T_ji / Z 
end

"""

    lianginfo_transfer(X::Matrix, dXdt::Matrix, dt::Real)

Compute the information transfer and its normalized value as defined by Liang (2021), as well as the Pearson's correlation.
"""
function lianginfo_transfer(X::Matrix, dXdt::Matrix, dt::Real)
    nvar, nt = size(X)

    # Compute covariance matrices and related quantities
    C = cov(X', X')
    dC = cov(X', dXdt')
    detC = det(C)
    invC = inv(C)
    Delta = Matrix(invC' .* detC)

    # Compute information transfer and Pearson's correlation
    T = zeros(nvar, nvar)
    R = copy(T)
    @inbounds for i in 1:nvar, j in 1:nvar
        T[j, i] = liang_index(detC, Delta[j,:], dC[:,i], C[i,j], C[i,i])
        R[j, i] = C[i,j] / sqrt(C[i,i] * C[j,j])
    end

    # Compute noise terms for normalized information transfer
    g = zeros(nvar)
    @inbounds for i in eachindex(g)
        a1k = invC * dC[:,i]
        f1 = mean(dXdt[i,:])
        @inbounds for k in eachindex(a1k)
            f1 -= a1k[k] * mean(X[k,:])
        end
        R1 = dXdt[i,:] .- f1
        @inbounds for k in eachindex(a1k)
            R1 .-= a1k[k] .* X[k, :]
        end
        Q1 = sum(R1 .^ 2)
        g[i] = Q1 * dt / nt
    end

    # Compute normalized information transfer
    tau = zeros(nvar, nvar)
    @inbounds for i in 1:nvar, j in 1:nvar
        tau[j,i] = liang_normindex(detC, Delta[i,:], dC[:,i], T[:,i], T[i,i], g[i], C[i,i], T[j,i])
    end
    return T, tau, R
end

"""

    bootstrapped_lianginfo_transfer(X::Matrix, dt::Real,
n_bootstrap::Int; k_euler::Int = 1)

Perform bootstrapping over [lianginfo_transfer](@ref).
"""
function bootstrapped_lianginfo_transfer(X::Matrix, t::AbstractVector,
    n_bootstrap::Int; k_euler::Int = 1)

    nvar, nt = size(X)
    dt = mean(diff(t))

    # Compute time derivative of states
    dXdt = forward_euler(X, t, k_euler)

    # Compute transfer information for original time series
    T, tau, R = lianginfo_transfer(X, dXdt, dt)

    # Initialize matrices for bootstrapping
    T_bootstrap = zeros(size(T)..., n_bootstrap)
    tau_bootstrap = zeros(size(tau)..., n_bootstrap)
    R_bootstrap = zeros(size(R)..., n_bootstrap)

    # Perform bootstrapping
    sgen = [surrogenerator(X[i, :], RandomFourier(), Xoshiro(123)) for i in 1:nvar]
    @threads for l in 1:n_bootstrap
        X_bootstrap = Matrix(hcat([sgen[i]() for i in 1:nvar]...)')
        dXdt_bootstrap = forward_euler(X_bootstrap, t, k_euler)
        T_, tau_, R_ = lianginfo_transfer(X_bootstrap, dXdt_bootstrap, dt)
        T_bootstrap[:, :, l] .= T_
        tau_bootstrap[:, :, l] .= tau_
        R_bootstrap[:, :, l] .= R_
    end

    # Compute sampled standard-deviation
    error_T = std(T_bootstrap, dims = 3)[:, :, 1]
    error_tau = std(tau_bootstrap, dims = 3)[:, :, 1]
    error_R = std(R_bootstrap, dims = 3)[:, :, 1]

    return InformationTransfer(T, tau, R, error_T, error_tau, error_R)
end