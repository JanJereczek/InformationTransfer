include("intro.jl")

function grace_eof(; interpolate = true)
    x, y, t, lat, lon, massbalance = load_aisgmb()
    nx, ny, nt = size(massbalance)

    if interpolate
        start_year = floor(t[1])
        start_month = ceil(get_month_from_decimaltime(t[1]))
        end_year = floor(t[end])
        end_month = floor(get_month_from_decimaltime(t[end]))
        t_even = collect(range( start_year + start_month/12, step = 1/12,
            stop = end_year + end_month/12 ))

        vecmassbalance = [massbalance[:, :, k] for k in axes(massbalance, 3)]
        itp = linear_interpolation(t, vecmassbalance)
        vecmassbalance_even = itp.(t_even)
        massbalance_even = cat(vecmassbalance_even..., dims = 3)
        anim_aisgmb(massbalance_even, filename = "AIS_GMB_even")
    end

    X, R, indices = reshape_without_missings(massbalance_even, t_even)
    Y, tcrop = moving_average(R, t_even, 6)
    U, s, V = tsvd(Y * Y', 100)
    plot_realtive_pc(s, "pc_ais_smoothed")
    return EOFresults(X, R, Y, "AIS mass balance", "kg/m2", s, U, t, indices, x, y, "none")
end

res = grace_eof()
pc1 = vec(F.U[:, 1]' * R)
pc2 = vec(F.U[:, 2]' * R)

eof1 = F.U[:, 1] .* F.Vt[:, 1]

eof1 = reshape_from_vec(indices, F.U[:, 1], nx, ny)
eof2 = reshape_from_vec(indices, F.U[:, 2], nx, ny)

function load_index(file::String)
    if (file !== "sam_index") && (file !== "enso34_index")
        error("Only accept sam_index and enso34_index as input.")
    end
    x_full = readdlm(datadir("exp_raw/$file.csv"), ',')
    x_vec = vcat([x_full[i, 2:end] for i in axes(x_full, 1)]...)
    yr = x_full[:, 1]
    t = vec(["$i-$j" for j in 1:12, i in yr])
    t = range(yr[1], step = 1/12, stop = yr[end]+11/12)
    return t, x_vec
end

# sam_full = readdlm(datadir("exp_raw/sam_index.csv"), ',')
# enso_full = readdlm(datadir("exp_raw/enso34_index.csv"), ',')

# sam_vec = vcat([sam_full[i, 2:end] for i in axes(sam_full, 1)]...)
# enso_vec = vcat([enso_full[i, 2:end] for i in axes(enso_full, 1)]...)

# yr_sam = sam_full[:, 1]
# t_sam_full = yr_sam[1]:1/12:yr_sam[end]+11/12

# yr_enso = enso_full[:, 1]
# t_enso_full = yr_enso[1]:1/12:yr_enso[end]+11/12

t_sam_full, sam_vec = load_index("sam_index")
t_enso_full, enso_vec = load_index("enso34_index")

i1_sam = findfirst(t_sam_full .>= first(t))-1
i1_enso = findfirst(t_enso_full .>= first(t))-1

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