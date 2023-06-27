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

agregate_oras(variablename = "sohtc700")
anim_oras(variablename = "sohtc700", stride = 1)