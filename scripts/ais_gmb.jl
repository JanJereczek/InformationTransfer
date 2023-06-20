include("intro.jl")
using NCDatasets, CairoMakie, TSVD

function reshape_without_missings(X3D::Array, t::Vector)
    nx, ny, nt = size(X3D)
    notmissing = sum(ismissing.(X3D), dims=3)[:, :, 1] .== 0
    indices = CartesianIndices(notmissing)[notmissing]
    X = zeros(length(indices), nt)
    R = zeros(length(indices), nt)

    for k in eachindex(indices)
        i, j = Tuple(indices[k])
        X[k, :] .= X3D[i, j, :]
        R[k, :] .= X3D[i, j, :] - linregression(X3D[i, j, :], t)
    end
    return X, R, indices
end

function linregression(x::Vector, t::Vector; lambda = 0.0)
    TT = hcat(t, ones(length(t)))'
    M = inv(TT * TT' + lambda .* LinearAlgebra.I(2) ) * TT
    m, p = M*x
    return m .* t .+ p
end

function pca(X::Matrix)
    C = cov(X)
end

function load_aisgmb()
    # https://data1.geo.tu-dresden.de/ais_gmb/
    ds = Dataset(datadir("exp_raw/AIS_GMB_grid.nc"))
    x = copy(ds["x"][:])
    y = copy(ds["y"][:])
    t = copy(ds["time_dec"][:])
    lat = copy(ds["lat"][:])
    lon = copy(ds["lon"][:])
    massbalance = copy(ds["dm"][:])
    close(ds)
    return x, y, t, lat, lon, massbalance
end

function anim_aisgmb(massbalance)
    nt = size(massbalance, 3)
    k = Observable(1)
    plot_mbalance = @lift(massbalance[:, :, $k])
    fig, ax, hm = heatmap(plot_mbalance, colorrange = (-1000, 1000), colormap = :balance)
    record(fig, "plots/AIS_GMB.mp4", 1:nt) do i
        k[] = i
    end
end

function reshape_eof(indices::Vector, EOF::Vector, nx::Int, ny::Int)
    X2D = zeros(nx, ny)
    for k in eachindex(indices)
        i, j = Tuple(indices[k])
        X2D[i, j] = EOF[k]
    end
    return X2D
end

x, y, t, lat, lon, massbalance = load_aisgmb()
nx, ny, nt = size(massbalance)
anim_aisgmb(massbalance)
X, R, indices = reshape_without_missings(massbalance, t)

# C = cov(R')
# U, S, V = svd(C)
# U * diagm(S) * V' ≈ C

F = svd(R * R')
# S[3] / sum(S)
pc1 = F.U[:, 1]' * R
pc2 = F.U[:, 2]' * R

eof1 = F.U[:, 1] .* F.Vt[:, 1]

eof1 = reshape_eof(indices, F.U[:, 1], nx, ny)
eof2 = reshape_eof(indices, F.U[:, 2], nx, ny)
heatmap(eof1)
heatmap(eof2)