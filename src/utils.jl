
function replace_missing!(X::Matrix{Union{Missing, T}}, x::Real) where {T <: AbstractFloat}
    x = T(x)
    mask = ismissing.(X)
    X[mask] .= x
end

# function replace_missing!(X::Matrix{Union{Missing, T}}) where {T <: AbstractFloat}
#     replace_missing_by_x!(X, 0)
# end

function replace_missing(X::Matrix{Union{Missing, T}}, x::Real) where {T <: AbstractFloat}
    Y  = copy(X)
    replace_missing!(Y, x)
    return Y
end

# function reshape_without_missings(X3D::Array, t::Vector)
#     nx, ny, nt = size(X3D)
#     notmissing = sum(ismissing.(X3D), dims=3)[:, :, 1] .== 0
#     indices = CartesianIndices(notmissing)[notmissing]
#     X = zeros(length(indices), nt)
#     R = zeros(length(indices), nt)

#     for k in eachindex(indices)
#         i, j = Tuple(indices[k])
#         X[k, :] .= X3D[i, j, :]
#         R[k, :] .= X3D[i, j, :] - linregression(X3D[i, j, :], t)
#     end
#     return X, R, indices
# end

function reshape_without_missings(X3D::Array, t::Vector)
    nx, ny, nt = size(X3D)
    mask = fill(true, nx, ny)
    return reshape_without_missings(X3D, t, mask)
end

function reshape_without_missings(X3D::Array, t::Vector, mask::BitMatrix)
    nx, ny, nt = size(X3D)
    notmissing = sum(ismissing.(X3D), dims=3)[:, :, 1] .== 0
    valid = mask .& notmissing
    indices = CartesianIndices(valid)[valid]
    X = zeros(length(indices), nt)
    R = zeros(length(indices), nt)

    for k in eachindex(indices)
        i, j = Tuple(indices[k])
        X[k, :] .= X3D[i, j, :]
        R[k, :] .= X3D[i, j, :] - linregression(X3D[i, j, :], t)
    end
    return X, R, indices
end


function linregression(x::Vector, t::AbstractVector; lambda = 0.0)
    TT = hcat(t, ones(length(t)))'
    M = inv(TT * TT' + lambda .* LinearAlgebra.I(2) ) * TT
    m, p = M*x
    return m .* t .+ p
end

function pca(X::Matrix)
    C = cov(X)
end

function plot_realtive_pc(s::Vector, plotname::String)
    relative_pc = s ./ sum(s)
    cumrelative_pc = cumsum(relative_pc)
    fig, ax, l = lines(relative_pc[1:100])
    lines!(ax, cumrelative_pc[1:100])
    save(plotsdir("$plotname.png"), fig)
    save(plotsdir("$plotname.pdf"), fig)
end

function reshape_eof(indices::Vector, EOF::Vector, nx::Int, ny::Int)
    X2D = Matrix{Union{Float64, Missing}}(fill(missing, nx, ny))
    for k in eachindex(indices)
        i, j = Tuple(indices[k])
        X2D[i, j] = EOF[k]
    end
    return X2D
end