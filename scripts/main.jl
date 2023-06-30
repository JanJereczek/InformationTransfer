include("intro.jl")

function get_principalcomponent(file)
    @load "$file" results
    return results.time, results.U' * results.Y
end

# function plot_eof_spatiotemporal(variablename::String; n_eof::Int = 2)
eof_dir = datadir("exp_pro/eof_oras")
files = readdir(eof_dir, join = true)
time, placeholder_pc = get_principalcomponent(files[1])
principal_components = zeros(eltype(placeholder_pc),
    size(placeholder_pc)..., length(files))
for k in axes(principal_components, 3)
    _, principal_components[:, :, k] = get_principalcomponent(files[k])
end

x, y, t, lat, lon, massbalance = load_aisgmb()
