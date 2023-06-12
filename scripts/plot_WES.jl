include("intro.jl")
include(srcdir("plotting.jl"))
using DelimitedFiles
using CairoMakie

X, lbl = readdlm(datadir("exp_raw/imbie_west_antarctica_2021_Gt.csv"), ',', header = true)
yr = X[:, 1]            # year
mb = X[:, 2]            # Mass balance (Gt/yr)
mb_sigma = X[:, 3]      # Mass balance uncertainty (Gt/yr)
cmb = X[:, 4]           # Cumulative mass balance (Gt)
cmb_sigma = X[:, 5]     # Cumulative mass balance uncertainty (Gt)

sam_full = readdlm(datadir("exp_raw/sam_index.csv"), ',')
enso_full = readdlm(datadir("exp_raw/enso34_index.csv"), ',')

i1_sam = findfirst(sam_full[:, 1] .>= first(yr))
i1_enso = findfirst(enso_full[:, 1] .>= first(yr))
i2_sam = findfirst(sam_full[:, 1] .>= last(yr)) - 1
i2_enso = findfirst(enso_full[:, 1] .>= last(yr)) - 1

sam = sam_full[i1_sam:i2_sam, :]
enso = enso_full[i1_enso:i2_enso, :]
sam_vec = vcat([sam[i, 2:end] for i in axes(sam, 1)]...)[eachindex(yr)]
enso_vec = vcat([enso[i, 2:end] for i in axes(enso, 1)]...)[eachindex(yr)]

ft = 30
lw = 5
yr_ticks = (1995:5:2020)
fig = Figure(resolution = (1600, 900), fontsize = ft)
nrow, ncol = 2, 2

function k(i::Int, j::Int)
    return (i-1) * ncol + j
end

xlabels = [
    "",
    "",
    L"Year $\,$",
    L"Year $\,$",
]

ylabels = [
    L"WAIS mass balance $\mathrm{(Gt \, yr^{-1})}$",
    L"ENSO-3.4 index $\,$",
    L"WAIS cum. mass balance (Gt) $\,$",
    L"SAM index $\,$",
]

xticklabelsvisible = [
    false,
    false,
    true,
    true,
]

axs = [Axis(
    fig[i, j],
    xlabel = xlabels[k(i, j)],
    ylabel = ylabels[k(i, j)],
    xticks = yr_ticks,
    xticklabelsvisible = xticklabelsvisible[k(i,j)],
) for i in 1:nrow, j in 1:ncol]

band!(axs[1], yr, mb-mb_sigma, mb+mb_sigma)
lines!(axs[1], yr, mb, color = :gray10, linewidth = lw)

band!(axs[2], yr, cmb-cmb_sigma, cmb+cmb_sigma)
lines!(axs[2], yr, cmb, color = :gray10, linewidth = lw)

barplot!(axs[3], yr, enso_vec, color = enso_vec, colormap = :balance)
barplot!(axs[4], yr, sam_vec, color = sam_vec, colormap = :balance)

save(plotsdir("WES.png"), fig)
save(plotsdir("WES.pdf"), fig)


dt = mean( diff(yr) )
X = Matrix(hcat(mb, sam_vec, enso_vec)')
n_bootstrap = 1000
dXdt = forward_euler(X, 1, dt)

infotransfer = bootstrapped_lianginfo_transfer(X, dt, n_bootstrap)

opts = (fontsize = 30, thin_lw = 2, thick_lw = 5)
fig = Figure(resolution = (800, 800), fontsize = opts.fontsize)
ax = Axis(fig[1,1])
vars = [L"$\dot{m}$", L"$i_\mathrm{SAM}$", L"$i_\mathrm{ENSO}$"]
plot_infotransfer!(ax, infotransfer.tau, infotransfer.error_tau, vars, opts)

# Introudce accronym WES (WAIS-ENSO-SAM)
save(plotsdir("infotransfer_WES.png"), fig)
save(plotsdir("infotransfer_WES.pdf"), fig)