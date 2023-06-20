#=
sam_itp = linear_interpolation(t_sam, sam_vec)
enso_itp = linear_interpolation(t_enso, enso_vec)

t_indices = t[1:end-1]
sam = sam_itp.(t_indices)
enso = enso_itp.(t_indices)
# i1_sam = findfirst(sam_full[:, 1] .>= first(t))
# i1_enso = findfirst(enso_full[:, 1] .>= first(t))

# sam = sam_full[i1_sam:end, :]
# enso = enso_full[i1_enso:end, :]
# sam_vec = vcat([sam[i, 2:end] for i in axes(sam, 1)]...)
# enso_vec = vcat([enso[i, 2:end] for i in axes(enso, 1)]...)

sam_sigma = cumsum(sam)
enso_sigma = cumsum(enso)

sam_filt = sam_sigma - linregression(sam_sigma, t_indices)
enso_filt = enso_sigma - linregression(enso_sigma, t_indices)
=#


# C = cov(R')
# U, S, V = svd(C)
# U * diagm(S) * V' ≈ C