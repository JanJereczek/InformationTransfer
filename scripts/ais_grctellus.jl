# https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.1_V3
dds = Dataset(datadir("exp_raw/GRCTellus.JPL.200204_202303.GLO.RL06.1M.MSCNv03CRI.nc"))
lwe = copy(dds["lwe_thickness"][:])
close(dds)
k = Observable(1)
plot_mbalance = @lift(lwe[:, :, $k])
fig, ax, hm = heatmap(plot_mbalance, colorrange = (-2000, 2000), colormap = :balance)
record(fig, "GRCTellus.mp4", 1:219) do i
    k[] = i
end
