include("../scripts/intro.jl")
include("../src/plotting.jl")

M = 2 .* (rand(5, 5) .- 0.5)
S = 0.5 .* (rand(5, 5) .- 0.5)
opts = (fontsize = 30, thin_lw = 2, thick_lw = 5)
fig = Figure(resolution = (800, 800), fontsize = opts.fontsize)
ax = Axis(fig[1,1])
vars = [L"a $\,$", L"b $\,$", L"c $\,$", L"d $\,$", L"e $\,$"]
plot_infotransfer!(ax, M, S, vars, opts)
# Colorbar(fig[1, 2], hmap; label = "values", width = 15, ticksize = 15)
fig