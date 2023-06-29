include("intro.jl")

function main()
    timestrings_osi, X_osi, lat_osi, lon_osi = load_osisaf(50)
    lat_oras, lon_oras, timestrings_oras, X_oras,
        units, longname = load_oras("ileadfra", 250)

    masking = "wais"

    if masking == "oras"
        lonlims = extrema(lon_oras)
        latlims = extrema(lat_oras)
    elseif masking == "wais"
        lonlims = (-180, -60)
        latlims = (-90, -60)
    end

    osi_mask = lonlatmask(lonlims, latlims, lon_osi, lat_osi)
    oras_mask = lonlatmask(lonlims, latlims, lon_oras, lat_oras)

    X_osi = cat(X_osi..., dims = 3) ./ 100
    nt = min(size(X_oras, 3), size(X_osi, 3))

    # X0_oras = Float32.(replace_missing(X_oras, 0f0))
    # X0_osi = Float32.(replace_missing(X_osi, 0f0))

    # mean_osi = [mean(X0_osi[osi_mask, k]) for k in axes(X0_osi, 3)[1:nt]]
    # mean_oras = [mean(X0_oras[oras_mask, k]) for k in axes(X0_oras, 3)[1:nt]]

    # fig, ax, l = lines(mean_osi)
    # lines!(ax, mean_oras)
    # fig

    # meanbit_osi = [mean(X0_osi[osi_mask, k] .> 0) for k in axes(X0_osi, 3)[1:nt]]
    # meanbit_oras = [mean(X0_oras[oras_mask, k] .> 0) for k in axes(X0_oras, 3)[1:nt]]

    # fig, ax, l = lines(meanbit_osi)
    # lines!(ax, meanbit_oras)
    # fig

    fighm = Figure()
    ax1 = Axis(fighm[1, 1], aspect = DataAspect())
    xlims!(ax1, (40, 180))
    ylims!(ax1, (100, 300))
    ax1.yreversed = true

    ax2 = GeoAxis(fighm[1, 2], source = "+proj=longlat +datum=WGS84",
        dest="+proj=sterea +lat_0=-90", aspect = AxisAspect(1/1.5))

    xlims!(ax2, lonlims...)
    ylims!(ax2, latlims...)

    for k in 1:nt
        println(k)

        contourf!(ax1, X_osi[:, :, k], colorrange = (0, 1), levels = 20)

        surface!(ax2, lon_oras[:, 1], lat_oras[1, :], X_oras[:, :, k],
            colorrange = (0, 1), shading = false)

        save(plotsdir("oras_vs_osisaf/$k.png"), fighm)
    end
end

main()