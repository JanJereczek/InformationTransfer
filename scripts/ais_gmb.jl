include("intro.jl")
include(srcdir("plotting.jl"))
using NCDatasets, CairoMakie, TSVD, DelimitedFiles, Interpolations

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

function linregression(x::Vector, t::AbstractVector; lambda = 0.0)
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
    X2D = Matrix{Union{Float64, Missing}}(fill(missing, nx, ny))
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

F = svd(R * R')
relative_pc = F.S ./ sum(F.S)
cumrelative_pc = cumsum(relative_pc)
fig, ax, l = lines(relative_pc[1:100])
lines!(ax, cumrelative_pc[1:100])
save(plotsdir("wes_pc.png"), fig)
save(plotsdir("wes_pc.pdf"), fig)

pc1 = vec(F.U[:, 1]' * R)
pc2 = vec(F.U[:, 2]' * R)

eof1 = F.U[:, 1] .* F.Vt[:, 1]

eof1 = reshape_eof(indices, F.U[:, 1], nx, ny)
eof2 = reshape_eof(indices, F.U[:, 2], nx, ny)

sam_full = readdlm(datadir("exp_raw/sam_index.csv"), ',')
enso_full = readdlm(datadir("exp_raw/enso34_index.csv"), ',')

sam_vec = vcat([sam_full[i, 2:end] for i in axes(sam_full, 1)]...)
enso_vec = vcat([enso_full[i, 2:end] for i in axes(enso_full, 1)]...)

yr_sam = sam_full[:, 1]
t_sam_full = yr_sam[1]:1/12:yr_sam[end]+11/12

yr_enso = enso_full[:, 1]
t_enso_full = yr_enso[1]:1/12:yr_enso[end]+11/12

i1_sam = findfirst(t_sam .>= first(t))-1
i1_enso = findfirst(t_enso .>= first(t))-1

sam_trunc = sam_vec[i1_sam:end]
enso_trunc = enso_vec[i1_enso:end]
t_sam = collect(t_sam_full[i1_sam:end])
t_enso = collect(t_enso_full[i1_enso:end])

sam_sigma = cumsum(sam_trunc)
enso_sigma = cumsum(enso_trunc)

sam_filt = sam_sigma - linregression(sam_sigma, t_sam)
enso_filt = enso_sigma - linregression(enso_sigma, t_enso)

sam_itp = linear_interpolation(t_sam, sam_filt, extrapolation_bc = Flat())
enso_itp = linear_interpolation(t_enso, enso_filt, extrapolation_bc = Flat())

normed(y::Vector) = y ./ maximum(abs.(y))

nrows, ncols = 2, 2
fig = Figure(resolution = (1600, 900), fontsize = 30)
axs = [Axis(fig[i, j]) for i in 1:nrows, j in 1:ncols]
lines!(axs[1], t_sam, normed(sam_filt), label = L"$\mathrm{SAM}_{\Sigma}$")
lines!(axs[1], t, normed(sam_itp.(t)), label = L"interpolated $\mathrm{SAM}_{\Sigma}$")
lines!(axs[1], t, -normed(pc1), label = L"Grace PC1 $\,$")
lines!(axs[2], t_enso, normed(enso_filt), label = L"$\mathrm{ENSO3.4}_{\Sigma}$")
lines!(axs[2], t, normed(enso_itp.(t)), label = L"interpolated $\mathrm{ENSO3.4}_{\Sigma}$")
lines!(axs[2], t, -normed(pc2), label = L"Grace PC2 $\,$")
xlims!(axs[1], (2002, 2022))
xlims!(axs[2], (2000, 2023))

crange = (-0.2, 0.2)
cmap = cgrad(:balance, rev = true)
heatmap!(axs[3], -eof1, colormap = cmap, colorrange = crange)
heatmap!(axs[4], -eof2, colormap = cmap, colorrange = crange)
axs[3].aspect = DataAspect()
axs[4].aspect = DataAspect()

fig
save(plotsdir("wes_eof.png"), fig)
save(plotsdir("wes_eof.pdf"), fig)

dt = mean( diff(t) )
X = Matrix(hcat(-normed(pc1), -normed(pc2), normed(sam_itp.(t)), normed(enso_itp.(t)))')
n_bootstrap = 1000
dXdt = forward_euler(X, t, 1)

infotransfer = bootstrapped_lianginfo_transfer(X, t, n_bootstrap)

opts = (fontsize = 30, thin_lw = 2, thick_lw = 5)
fig = Figure(resolution = (800, 800), fontsize = opts.fontsize)
ax = Axis(fig[1,1])
vars = [L"$\mathrm{PC1}$", L"$\mathrm{PC2}$",
    L"$\mathrm{SAM}_{\Sigma}$", L"$\mathrm{ENSO}_{\Sigma}$"]
plot_infotransfer!(ax, infotransfer.tau, infotransfer.error_tau, vars, opts)

# Introudce accronym WES (WAIS-ENSO-SAM)
save(plotsdir("infotransfer_wes_eof.png"), fig)
save(plotsdir("infotransfer_wes_eof.pdf"), fig)