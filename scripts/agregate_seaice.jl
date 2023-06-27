include("intro.jl")

function latlon_seaice(; k=50)
    
    testfile1 = datadir("osisaf.met.no/reprocessed/ice/conc/v2p0/2001/01/" *
        "ice_conc_sh_ease2-250_cdr-v2p0_200101011200.nc")
    ds = NCDataset(testfile1)
    lat = copy(ds["lat"][:][k:end-k+1, k:end-k+1])
    lon = copy(ds["lon"][:][k:end-k+1, k:end-k+1])
    close(ds)
    jldsave(datadir("exp_pro/osisaf/latlon_k=$(k).jld2"),
        lat = lat, lon = lon)
end

function agregate_seaice_monthly(; k = 50)
    testfile1 = datadir("osisaf.met.no/reprocessed/ice/conc/v2p0/1992/01/"
        * "ice_conc_sh_ease2-250_cdr-v2p0_199201011200.nc")
    ds = NCDataset(testfile1)
    lat1 = copy(ds["lat"][:][k:end-k+1, k:end-k+1])
    lon1 = copy(ds["lon"][:][k:end-k+1, k:end-k+1])
    ice_conc1 = copy(ds["ice_conc"][:][k:end-k+1, k:end-k+1])
    close(ds)

    testfile2 = datadir("osisaf.met.no/reprocessed/ice/conc-cont-reproc/v2p0/2016/01/"
        * "ice_conc_sh_ease2-250_icdr-v2p0_201601011200.nc")
    ds = NCDataset(testfile2)
    lat2 = copy(ds["lat"][:][k:end-k+1, k:end-k+1])
    lon2 = copy(ds["lon"][:][k:end-k+1, k:end-k+1])
    ice_conc2 = copy(ds["ice_conc"][:][k:end-k+1, k:end-k+1])
    close(ds)

    if (lat2 ≈ lat1) & (lon2 ≈ lon1) & (size(ice_conc1) == size(ice_conc2))
        println("Coherent grid across data sets.")
    else
        error("Incoherent grid across data sets.")
    end

    years = 1992:2022
    months = 1:12
    timestrings = vec(["$year-$month" for month in months, year in years])
    replace_missings = false
    if replace_missings
        month_averaged_seaice = [zeros(Float32, size(ice_conc1)...) for t in timestrings]
    else
        month_averaged_seaice = [zeros(Union{Missing, Float32}, size(ice_conc1)...) for t in timestrings]
    end

    dirs1992_2015 = readdir(datadir("osisaf.met.no/reprocessed/ice/conc/v2p0"), join = true)
    dirs2016_2022 = readdir(datadir("osisaf.met.no/reprocessed/ice/conc-cont-reproc/v2p0"), join = true)
    yeardirs = vcat(dirs1992_2015, dirs2016_2022)

    for i_year in eachindex(yeardirs)
        monthdirs = readdir(yeardirs[i_year], join=true)
        popfirst!(monthdirs)
        # println(size(monthdirs))

        for i_month in eachindex(monthdirs)
            daydirs = readdir(monthdirs[i_month], join=true)
            popfirst!(daydirs)
            j = (i_year-1)*12 + i_month
            println("$j:  $(length(daydirs))")  #  $(last(daydirs))")
            for i_day in eachindex(daydirs)
                NCDataset(daydirs[i_day]) do ds
                    if replace_missings
                        month_averaged_seaice[j] .+= replace_missing_by_zero!(
                            copy(ds["ice_conc"][:][k:end-k+1, k:end-k+1]))
                    else
                        month_averaged_seaice[j] .+= (copy(
                            ds["ice_conc"][:][k:end-k+1, k:end-k+1]))
                    end
                end
            end
            month_averaged_seaice[j] ./= length(daydirs)
        end
    end
    jldsave(datadir("exp_pro/osisaf/monthly_agregated_k=$(k).jld2"),
        timestrings = timestrings,
        month_averaged_seaice = month_averaged_seaice)
end

function anim_seaice(; k = 50)
    data = jldopen(datadir("exp_pro/osisaf/monthly_agregated_k=$(k).jld2"))
    timestrings, month_averaged_seaice = data["timestrings"], data["month_averaged_seaice"]

    j = Observable(1)
    observable_seaice = @lift( month_averaged_seaice[$j] )
    observable_timestring = @lift( "Time: $(timestrings[$j])" )

    ft = 30
    cmap = cgrad(:Blues, rev = true)
    fig = Figure(resolution = (1300, 1200), fontsize = ft)
    ax = Axis(fig[1,1], aspect = DataAspect(), yreversed = true)
    hidedecorations!(ax)

    crange = (0, 100)
    heatmap!(ax, observable_seaice, colormap = cmap, colorrange = crange)
    text!(ax, 150, 150, text = observable_timestring, fontsize = ft)
    Colorbar(fig[1, 2], height = Relative(0.6), colorrange = crange, colormap = cmap,
        label = L"Sea-ice fraction (%) $\,$")
    contour!(ax, replace_missing(month_averaged_seaice[1], -1), levels = [-0.5])
    record(fig, plotsdir("osisaf/monthly_averaged_k=$k.mp4"), eachindex(timestrings), framerate = 10) do jj
        j[] = jj
    end
end

# agregate_seaice_monthly()
# anim_seaice()


# function mask_seaice(; k = 50)
k = 50
data = jldopen(datadir("exp_pro/osisaf/monthly_agregated_k=$(k).jld2"))
timestrings = data["timestrings"]
month_averaged_seaice = data["month_averaged_seaice"]
latlon = jldopen(datadir("exp_pro/osisaf/latlon_k=$(k).jld2"))
lat, lon = latlon["lat"], latlon["lon"]
latlims = (-80, -60)
lonlims = (-180, -60)
mask = (latlims[1] .< lat .< latlims[2]) .& (lonlims[1] .< lon .< lonlims[2])
fig, ax, hm = heatmap(month_averaged_seaice[1])
contour!(ax, mask, levels = [0.5])
fig

X3D = cat(month_averaged_seaice..., dims=3)
X, R, indices = reshape_without_missings(X3D, collect(eachindex(timestrings)), mask)
# TODO: subtract seasonnality

U, s, V = tsvd(R * R', 100)
plot_realtive_pc(s, "osisaf/pc")

pc1 = vec(U[:, 1]' * R)
# eof1 = U[:, 1] .* V[:, 1]'
nx, ny, nt = size(X3D)
eof1 = reshape_eof(indices, U[:, 1], nx, ny)
heatmap(eof1)