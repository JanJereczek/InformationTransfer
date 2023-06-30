include("intro.jl")

x, y, t, lat, lon, massbalance = load_aisgmb()
nx, ny, nt = size(massbalance)
X, R, indices = reshape_without_missings(massbalance, t)

F = svd(R * R')

plot_realtive_pc(F.S, "grace/wes_svd_S")

pcs = F.U[:, 1:2]' * R
t = collect(axes(pcs, 2))
dpcsdt = forward_euler(pcs, t, 1)
smooth_pcs, t_smooth = moving_average(pcs, t, 6)
smooth_dpcsdt = forward_euler(smooth_pcs, t_smooth, 1)

fig = Figure()
axs = [Axis(fig[i, 1]) for i in 1:2]

for i in 1:2
    lines!(axs[1], t, pcs[i, :], color = Cycled(i))
    lines!(axs[1], t_smooth, smooth_pcs[i, :], color = Cycled(i), linestyle = :dash)
    lines!(axs[2], t, dpcsdt[i, :], color = Cycled(i))
    lines!(axs[2], t_smooth, smooth_dpcsdt[i, :], color = Cycled(i), linestyle = :dash)
end
fig
