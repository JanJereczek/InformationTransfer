
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
