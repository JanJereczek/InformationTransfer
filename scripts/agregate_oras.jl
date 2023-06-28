include("intro.jl")

function agregate_oras(; variablename = "sohtc300", k = 250)
    variabledir = datadir("oras5/ORCA025/$(variablename)/opa0")
    yeardirs = readdir(variabledir, join = true)

    testfile = datadir("$variabledir/$(variablename)_ORAS5_1m_1992/" *
    "$(variablename)_ORAS5_1m_199201_grid_T_02.nc")
    ds = NCDataset(testfile)
    full_lat = copy(ds["nav_lat"][:])
    full_lon = copy(ds["nav_lon"][:])
    full_var = copy(ds[variablename][:])
    units = ds[variablename].attrib["units"]
    longname = ds[variablename].attrib["long_name"]
    close(ds)

    years = 1992:2018
    months = 1:12
    timestrings = vec(["$year-$month" for month in months, year in years])

    lat = full_lat[:, 1:k]
    lon = full_lon[:, 1:k]
    toolvar = full_var[:, 1:k]
    stacked_var = [zeros(eltype(toolvar), size(toolvar)...) for t in timestrings]

    for i_year in eachindex(yeardirs)
        files = readdir(yeardirs[i_year], join = true)
        for i_month in eachindex(files)
            j = (i_year-1)*12 + i_month
            NCDataset(files[i_month]) do ds
                stacked_var[j] .= replace_missing(ds[variablename][:][:, 1:k], -1f10)
            end
        end
    end
    jldsave(datadir("exp_pro/oras5/$(variablename)_k=$k.jld2"), lat = lat, lon = lon,
        timestrings = timestrings, stacked_var = stacked_var, units = units,
        longname = longname)
end

function anim_oras(; variablename = "sohtc300", k = 250, stride = 1)
    data = jldopen(datadir("exp_pro/oras5/$(variablename)_k=$k.jld2"))
    lat, lon, timestrings = data["lat"], data["lon"], data["timestrings"]
    var, units, longname = data["stacked_var"], data["units"], data["longname"]
    
    cmin = minimum([minimum(var[i]) for i in eachindex(var)])
    cmax = maximum([maximum(var[i]) for i in eachindex(var)])
    crange = (cmin, cmax)

    ftsize = 30
    cmap = cgrad(:jet)
    fig = Figure(resolution = (1300, 1200), font = srcdir("cmunrm.ttf"),
        fontsize = ftsize)
    ga = GeoAxis(fig[1, 1], source = "+proj=longlat +datum=WGS84",
        dest="+proj=sterea +lat_0=-90")        # coastlines = true,
    ga.aspect = AxisAspect(1)
    xlims!(ga, -180, 180)
    ylims!(ga, -90, maximum(lat))

    j = Observable(1)
    plotvar = @lift(var[$j])
    observable_timestring = @lift( "Time: $(timestrings[$j])" )

    text!(ga, 0.0, -90, text = observable_timestring, fontsize = ftsize)
    sf = surface!(ga, lon, lat, plotvar; shading = false, colormap = cmap,
        colorrange = crange, lowclip = :transparent)
    Colorbar(fig[1, 2], sf, height = Relative(0.5), label = L"%$longname (%$units) $\,$")
    # hidedecorations!(ga)
    record(fig, plotsdir("oras/$(variablename)_stride=$(stride)_k=$k.mp4"), 
        eachindex(timestrings)[1:stride:end], framerate = 10) do jj
        j[] = jj
    end
    # fig
end

function agregate_oras_full(; variablename = "sohtc300", k = 250)
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

function anim_oras_full(; variablename = "sohtc300", k = 250, stride = 1)
    data = jldopen(datadir("exp_pro/oras5/$(variablename)_k=$k.jld2"))
    lat, lon, timestrings = data["lat"], data["lon"], data["timestrings"]
    X, units, longname = data["X"], data["units"], data["longname"]

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

vnames = ["ileadfra", "sohtc300", "sohtc700", "sohtcbtm", "sometauy", "sosaline",
    "sowaflup", "sozotaux"]

agregate_data = false
plot_data = true

if agregate_data
    for vname in vnames[7:end]
        agregate_oras_full(variablename = vname)
    end
end

if plot_data
    vname = vnames[1]
    anim_oras_full(variablename = vname, stride = 2)
end