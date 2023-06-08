include("intro.jl")
using DelimitedFiles
using CairoMakie

X, lbl = readdlm(datadir("exp_raw/imbie_west_antarctica_2021_Gt.csv"), ',', header = true)
yr = X[:, 1]            # year
mb = X[:, 2]            # Mass balance (Gt/yr)
mb_sigma = X[:, 3]      # Mass balance uncertainty (Gt/yr)
cmb = X[:, 4]           # Cumulative mass balance (Gt)
cmb_sigma = X[:, 5]     # Cumulative mass balance uncertainty (Gt)

ft = 30
lw = 5
yr_ticks = (1995:5:2020)
fig = Figure(resolution = (1600, 900), fontsize = ft)
nrow, ncol = 2, 1

function k(i::Int, j::Int)
    return (i-1) * ncol + j
end

xlabels = [
    "",
    L"Year $\,$",
]

ylabels = [
    L"WAIS mass balance $\mathrm{(Gt \, yr^{-1})}$",
    L"WAIS cum. mass balance (Gt) $\,$",
]

xticklabelsvisible = [
    false,
    true,
]

axs = [Axis(
    fig[i, j],
    xlabel = xlabels[k(i,j)],
    ylabel = ylabels[k(i,j)],
    xticks = yr_ticks,
    xticklabelsvisible = xticklabelsvisible[k(i,j)],
) for i in 1:nrow, j in 1:ncol]

band!(axs[1], yr, mb-mb_sigma, mb+mb_sigma)
lines!(axs[1], yr, mb, color = :gray10, linewidth = lw)

band!(axs[2], yr, cmb-cmb_sigma, cmb+cmb_sigma)
lines!(axs[2], yr, cmb, color = :gray10, linewidth = lw)
fig