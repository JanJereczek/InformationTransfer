using CairoMakie

"""

    rectangle(i, j)

Return a `Point2f` containing the vertices of the bounding box for the text at
`(i, j)` of the heatmap.
"""
function rectangle(i::Int, j::Int)
    return Point2f[(i-0.4, j-0.4), (i-0.4, j+0.4), (i+0.4, j+0.4), (i+0.4, j-0.4)]
end

"""

    plot_infotransfer!(ax, M, S, vars, opts)

Plot onto `ax` the annotated heatmap of `M::Matrix` = a statistical metric.
Add bounding box if a value is significant given the standard-deviation `S::Matrix`.
Annotate axes with `vars` and set plotting options through `opts::NamedTuple`.
"""
function plot_infotransfer!(ax, M::Matrix, S::Matrix, vars::Vector, opts::NamedTuple;
    metric::String = "Normalized information transfer")

    n1, n2 = size(M)
    m1, m2 = size(S)
    if n1 != n2
        error("M must be a square matrix.")
    elseif (n1 != m1) || (n2 != m2)
        error("S must have same dimensions as M.")
    elseif length(vars) != n1
        error("Labels contained in vars must have same length as M.")
    end

    if metric == "Normalized information transfer"
        crange = (-100, 100)
    elseif metric == "Pearson's correlation coefficient"
        crange = (-1, 1)
    else
        error("The chosen metric must be either the normalized information" *
            "transfer or Pearson's correlation coefficient")
    end
    heatmap!(ax, M, colormap = :balance, colorrange = crange)

    for idx in CartesianIndices(M)
        i, j = Tuple(idx)
        txtcolor = abs(M[i, j]) < 25.0 ? :black : :white
        text!(ax, L"%$(round(M[i,j], digits = 2)) $\,$", position = (i, j), color = txtcolor,
            align = (:center, :center), fontsize = opts.fontsize)
        
        # Only draw rectangle around if significant
        if (M[i, j] + S[i, j]) * (M[i, j] - S[i, j]) > 0
            poly!(rectangle(i, j), color = :transparent,
                strokecolor = txtcolor, strokewidth = opts.thick_lw)
        end
    end
    ax.title = L"$\mathrm{\textbf{%$metric}...}$"
    ax.xlabel = L"$\mathrm{\textbf{To}...}$"
    ax.ylabel = L"$\mathrm{\textbf{From}...}$"

    ticks = ( 1:length(vars), vars )
    ax.xticks = ticks
    ax.yticks = ticks
    ax.xticklabelrotation = π / 4
    ax.xticklabelalign = (:right, :center)

    ax.yreversed = true
end

function anim_aisgmb(massbalance; filename::String = "AIS_GMB")
    nt = size(massbalance, 3)
    maxabs = maximum(abs.(skipmissing(massbalance)))
    crange = (-maxabs, maxabs)
    cmap = cgrad(:balance, rev = true)
    fig = Figure(resolution = (1300, 900), fontsize = 30)
    ax = Axis(fig[1, 1], aspect = DataAspect())
    hidedecorations!(ax)
    Colorbar(fig[1, 2], height = Relative(0.5), colormap = cmap,
        colorrange = crange, label = L"Mass balance $\mathrm{kg \, m^{-2}}$")

    k = Observable(1)
    plot_mbalance = @lift(massbalance[:, :, $k])
    heatmap!(ax, plot_mbalance, colorrange = crange, colormap = cmap, lowclip = :transparent)
    record(fig, "plots/$filename.mp4", 1:nt) do i
        k[] = i
    end
end