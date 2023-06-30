#######################################################
# Data manipulation
#######################################################

function replace_missing(X::Array{Union{Missing, T}}, x::Real) where {T<:AbstractFloat}
    Y = copy(X)
    replace_missing!(Y, x)
    return Y
end

function replace_missing!(X::Array{Union{Missing, T}}, x::Real) where {T<:AbstractFloat}
    x = T(x)
    mask = ismissing.(X)
    X[mask] .= x
end

function reshape_without_missings(X3D::Array, t::Vector; q::Real = 1)
    nx, ny, nt = size(X3D)
    mask = BitMatrix(fill(true, nx, ny))
    return reshape_without_missings(X3D, t, mask, q)
end

function reshape_without_missings(X3D::Array{Union{Missing, T}},
    t::Vector, mask::BitMatrix, q::Real) where {T<:AbstractFloat}
    nx, ny, nt = size(X3D)
    notmissing = sum(ismissing.(X3D), dims=3)[:, :, 1] .== 0
    valid = mask .& notmissing
    valid_indices = CartesianIndices(valid)[valid]
    keep = randomprune_indices(valid_indices, q)
    indices = valid_indices[keep]

    X = zeros(T, length(indices), nt)
    R = zeros(T, length(indices), nt)

    for k in eachindex(indices)
        i, j = Tuple(indices[k])
        X[k, :] .= X3D[i, j, :]
        R[k, :] .= X3D[i, j, :] - linregression(T.(X3D[i, j, :]), t)
    end
    return X, R, indices
end

function reshape_from_vec(indices::Vector, v::Vector, nx::Int, ny::Int)
    X2D = Matrix{Union{Float64, Missing}}(fill(missing, nx, ny))
    for k in eachindex(indices)
        i, j = Tuple(indices[k])
        X2D[i, j] = v[k]
    end
    return X2D
end

function recover_lonlat_from_indices(lon::Matrix, lat::Matrix, indices::Vector{CartesianIndex{2}})
    lonvec = zeros(eltype(lon), length(indices))
    latvec = zeros(eltype(lon), length(indices))
    @inbounds for idx in eachindex(indices)
        i, j = indices[idx].I
        lonvec[idx] = lon[i, j]
        latvec[idx] = lat[i, j]
    end
    return lonvec, latvec
end

function randomprune_indices(x::Vector, q::Real)
    return rand(length(x)) .<= q
end

get_month_from_decimaltime(t::AbstractFloat) = (t - floor(t)) * 12

#######################################################
# Filtering
#######################################################

function linregression(x::Vector{T}, t::AbstractVector; lambda = T(0)) where {T<:Real}
    TT = hcat(t, ones(T, length(t)))'
    M = inv(TT * TT' + lambda .* LinearAlgebra.I(2) ) * TT
    m, p = M*x
    return m .* t .+ p
end

function moving_average(X::Array{T, 3}, window_hw::Int) where {T<:Real}
    Y = zeros(T, size(X)...)
    for k in axes(X, 3)[window_hw+1:end-window_hw]
        Y[:, :, k] .= mean(X[:, :, k-window_hw:k+window_hw-1], dims=3)
    end
    return Y
end

function moving_average(X::Matrix{T}, t::AbstractVector, window_hw::Int) where {T<:Real}
    
    tcrop = t[window_hw+1:end-window_hw]
    Y = zeros(T, size(X,1), length(tcrop))
    for k in axes(Y, 2)
        Y[:, k] .= mean(X[:, k:k+2*window_hw-1], dims=2)
    end
    return Y, tcrop
end

#######################################################
# Plotting
#######################################################

function plot_realtive_pc(s::Vector, plotname::String)
    relative_pc = s ./ sum(s)
    cumrelative_pc = cumsum(relative_pc)

    nplot = findfirst(cumrelative_pc .> 0.99)
    ftsize = 30
    fig = Figure(resolution = (1600, 900), fontsize = ftsize)
    ax = Axis(fig[1, 1], xlabel = L"EOF number $\,$", ylabel = L"EOF weight $\,$")
    colors = [:gray, :cornflowerblue]
    barplot!(cumrelative_pc[1:nplot], color = colors[1], strokecolor = :black, strokewidth = 1)
    barplot!(relative_pc[1:nplot], color = colors[2], strokecolor = :black, strokewidth = 1)

    # Legend
    labels = [L"individual $\,$", L"cumulative $\,$"]
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    title = ""
    Legend(fig[1,2], elements, labels, title)

    # xlims!(ax, (-1, nmax+1))
    save(plotsdir("$plotname.png"), fig)
    save(plotsdir("$plotname.pdf"), fig)
end

#######################################################
# Data loading
#######################################################

function load_aisgmb()
    # https://data1.geo.tu-dresden.de/ais_gmb/
    ds = Dataset(datadir("exp_raw/AIS_GMB_grid.nc"))
    x = copy(ds["x"][:])
    y = copy(ds["y"][:])
    t = copy(ds["time_dec"][:])
    lat = copy(ds["lat"][:])
    lon = copy(ds["lon"][:])
    massbalance = copy(ds["dm"][:])
    close(ds)
    return x, y, t, lat, lon, massbalance
end

function load_oras(variablename::String, k::Int)
    data = jldopen(datadir("exp_pro/oras5/$(variablename)_k=$k.jld2"))
    lat, lon, timestrings = data["lat"], data["lon"], data["timestrings"]
    var, units, longname = data["X"], data["units"], data["longname"]
    return lat, lon, timestrings, var, units, longname
end

function load_osisaf(k)
    data = jldopen(datadir("exp_pro/osisaf/monthly_agregated_k=$(k).jld2"))
    latlon = jldopen(datadir("exp_pro/osisaf/latlon_k=$(k).jld2"))
    timestrings = data["timestrings"]
    month_averaged_seaice = data["month_averaged_seaice"]
    lat, lon = latlon["lat"], latlon["lon"]
    return timestrings, month_averaged_seaice, lat, lon
end

#######################################################
# Masking
#######################################################

function lonlatmask(lonlims, latlims, lon, lat)
    return (latlims[1] .< lat .< latlims[2]) .& (lonlims[1] .< lon .< lonlims[2])
end

acc_mask = (name = "ACC", latlims = (-62, -58), lonlims = (-140, -60))
amundsen_mask = (name = "Amundsen", latlims = (-90, -60), lonlims = (-140, -60))
polynia_mask = (name = "Polynia", latlims = (-80, -70), lonlims = (-140, -60))

#######################################################
# Structs
#######################################################
struct EOFresults{T<:AbstractFloat}
    X2D::Matrix{T}                      # original data (reshaped)
    R::Matrix{T}                        # detrended data (reshaped)
    Y::Matrix{T}                        # (annuly) smoothed data (reshaped)
    longname::String                    # for plotting
    units::String                       # for plotting
    s::Vector{T}                        # from (truncated) svd
    U::Matrix{T}                        # from (truncated) svd
    time::Vector                        # time, either numerical or string
    indices::Vector{CartesianIndex{2}}  # indices of cells that have no missing
    x::Vector{T}                        # either x from projection or longitude
    y::Vector{T}                        # either y from projection or latitude
    gridtype::String                    # keyword to know which projection
end