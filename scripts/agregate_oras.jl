include("intro.jl")

function agregate_oras(; variablename = "sohtc300", k = 250)
    variabledir = datadir("oras5_full/$(variablename)")

    testfile = "$variabledir/$(variablename)_control_monthly_highres_2D_199201_CONS_v0.1.nc"
    ds = NCDataset(testfile)
    full_lat = copy(ds["nav_lat"][:])
    full_lon = copy(ds["nav_lon"][:])
    full_var = copy(ds[variablename][:][:, :, 1])
    units = ds[variablename].attrib["units"]
    longname = ds[variablename].attrib["long_name"]
    close(ds)

    monthfiles = readdir(variabledir, join = true)
    nt = length(monthfiles)
    years = 1992:2023
    months = 1:12
    timestrings = vec(["$year-$month" for month in months, year in years])[1:nt]

    lat = Float32.(full_lat[:, 1:k])
    lon = Float32.(full_lon[:, 1:k])
    toolvar = full_var[:, 1:k]
    X = zeros(eltype(toolvar), size(toolvar)..., nt)

    for i_month in eachindex(monthfiles)
        NCDataset(monthfiles[i_month]) do ds
            println("Agregating $variablename for t=$(timestrings[i_month])")
            X[:, :, i_month] .= copy(ds[variablename][:][:, 1:k])
        end
    end
    jldsave(datadir("exp_pro/oras5/$(variablename)_k=$k.jld2"), lat = lat, lon = lon,
        timestrings = timestrings, X = X, units = units, longname = longname)
end

function anim_oras(; variablename = "sohtc300", k = 250, stride = 1)
    lat, lon, timestrings, X, units, longname = load_oras(variablename, k)

    latlims = (-90, -60)
    lonlims = (-180, -60)
    latrange, lonrange = latlims[2]-latlims[1], lonlims[2]-lonlims[1]

    t = collect(eachindex(timestrings))
    mask = (latlims[1] .< lat .< latlims[2]) .& (lonlims[1] .< lon .< lonlims[2])
    Y, R, indices = reshape_without_missings(X, t, mask)
    crange = extrema(Y)

    latvec, lonvec = lat[1, :], lon[:, 1]
    latmask = latlims[1] .< latvec .< latlims[2]
    lonmask = lonlims[1] .< lonvec .< lonlims[2]
    Z = X[lonmask, latmask, :]

    ftsize = 30
    cmap = cgrad(:jet)
    fig = Figure(resolution = (900, 1200), font = srcdir("cmunrm.ttf"),
        fontsize = ftsize)
    ga = GeoAxis(fig[1, 1], source = "+proj=longlat +datum=WGS84",
        dest="+proj=sterea +lat_0=-90")        # coastlines = true,
    ga.aspect = AxisAspect(1/(1+sin(deg2rad(lonrange-90))))# lonrange/latrange)
    xlims!(ga, lonlims...)
    ylims!(ga, latlims...)

    j = Observable(1)
    plotvar = @lift(Z[:, :, $j])
    observable_timestring = @lift( "Time: $(timestrings[$j])" )

    text!(ga, -100.0, -80, text = observable_timestring, fontsize = ftsize)
    cf = contourf!(ga, lonvec[lonmask], latvec[latmask], plotvar;
        shading = false, colormap = cmap,  lowclip = :transparent,
        levels = range(crange[1], stop = crange[end], length = 50))
    Colorbar(fig[1, 2], cf, height = Relative(0.5), label = L"%$longname (%$units) $\,$")

    println("Animating $variablename...")
    record(fig, plotsdir("oras/$(variablename)_stride=$(stride)_k=$k.mp4"), 
        eachindex(timestrings)[1:stride:end], framerate = 10) do jj
        println("Frame $jj")
        j[] = jj
    end
    println("Finished!")
    println("---------------------")
end

# Pass mask as argument: ACC, Amundsen, Polynia, WAISKatabatic
function eof_oras_ocean(variablename; q = 0.25, nsv = 100)
    k = 250
    maskname = "Amundsen"

    println("------------------------")
    println("$variablename:")
    println("   Loading data...")
    lat, lon, timestrings, X3D, units, longname = load_oras(variablename, k)
    nx, ny, nt = size(X3D)

    println("   Preprocessing data...")
    latlims = (-90, -60)
    lonlims = (-140, -60)
    latrange, lonrange = latlims[2]-latlims[1], lonlims[2]-lonlims[1]

    t = collect(eachindex(timestrings))
    mask = lonlatmask(lonlims, latlims, lon, lat)
    X2D, R, indices = reshape_without_missings(X3D, t, mask, q)
    lonvec, latvec = recover_lonlat_from_indices(lon, lat, indices)

    Y, tcrop = moving_average(R, timestrings, 6)
    println("   Computing SVD...")
    U, s, V = tsvd(Y * Y', nsv)

    println("   Saving results...")
    results = EOFresults(X2D, R, Y, longname, units, s, U, tcrop, indices,
        lonvec, latvec, "lonlatvec")
    jldsave(datadir("exp_pro/eof_oras/$(variablename)_$(maskname).jld2"),
        results = results)
end

function plot_eof_weights(variablename::String)
    maskname = "Amundsen"
    file = datadir("exp_pro/eof_oras/$(variablename)_$(maskname).jld2")
    @load "$file" results
    plot_realtive_pc(results.s, "oras/pc_$(variablename)_$(maskname)")
end

function plot_eof_spatiotemporal(variablename::String; n_eof::Int = 2)
    maskname = "Amundsen"
    file = datadir("exp_pro/eof_oras/$(variablename)_$(maskname).jld2")
    @load "$file" results
    longname, units, s, U = results.longname, results.units, results.s, results.U
    lonvec, latvec = results.x .+ 100, results.y
    t = results.time
    
    nsv = length(s)
    pcs = U' * results.Y
    
    ftsize = 60
    margin = 0.25
    msize = 15
    lw = 10
    xlims = extrema(lonvec) .+ (-margin, margin)
    ylims = extrema(latvec) .+ (-margin, margin)
    crange = extrema(U)
    cmap = cgrad(:jet)
    xtks_stride = 24
    xtks = (eachindex(t)[1:xtks_stride:end], t[1:xtks_stride:end])

    aspectratio = 0.35*(xlims[2]-xlims[1])/(ylims[2]-ylims[1])

    # println(typeof(lonvec), typeof(latvec), typeof(U[:, 1]))
    fig = Figure(resolution = (3200, 2000), font = srcdir("cmunrm.ttf"),
        fontsize = ftsize)
    Colorbar(fig[1, 1:n_eof], colormap = cmap, colorrange = crange,
        vertical = false, width = Relative(0.3), height = Relative(0.1), label = L"%$longname (1)$\,$")
    for j in 1:n_eof
        ga = GeoAxis(fig[2:4, j], source = "+proj=longlat +datum=WGS84",
            dest="+proj=sterea +lat_0=-90", title = L"EOF n°%$j $\,$",
            aspect = AxisAspect(aspectratio))
        scatter!(ga, lonvec, latvec; color = U[:, j], 
            shading = false, colormap = cmap, markersize = msize,
            levels = range(crange[1], stop = crange[end], length = 50))
        xlims!(ga, xlims)
        ylims!(ga, ylims)
        ax = Axis(fig[5, j], xlabel = L"Time $\,$", xticks = xtks,
            xticklabelrotation = π/2)
        lines!(ax, pcs[j, :], linewidth = lw)
        if j == 1
            ax.ylabel = L"EOF magnitude (%$units) $\,$"
        end
    end
    rowgap!(fig.layout, 1, -200)
    figname = "eof_$variablename"
    save(plotsdir("oras/$figname.pdf"), fig)
    save(plotsdir("oras/$figname.png"), fig)
end

# function eof_oras_wind()
# end

global VARIABLENAMES = ["ileadfra", "sohtc300", "sohtc700", "sohtcbtm",
    "sosaline", "sowaflup", "sozotaux", "sometauy"]

global COLORMAPS = [ cgrad(:Blues, rev = true), :jet, :jet, :jet,
    :balance, :thermal, missing, missing ]

agregate_data_bool = false
plot_data_bool = false
compute_eof_bool = true
plot_eof_spatiotemporal_bool = true
plot_eof_weights_bool = true
for variable in VARIABLENAMES[end:end]
    agregate_data_bool ? agregate_oras(variablename = variable) : nothing
    plot_data_bool ? anim_oras(variablename = variable, stride = 3) : nothing
    compute_eof_bool ? eof_oras_ocean(variable) : nothing
    plot_eof_spatiotemporal_bool ? plot_eof_spatiotemporal(variable) : nothing
    plot_eof_weights_bool ? plot_eof_weights(variable) : nothing
end