using CairoMakie

function rectangle(i::Int, j::Int)
    return Point2f[(i-0.4, j-0.4), (i-0.4, j+0.4), (i+0.4, j+0.4), (i+0.4, j-0.4)]
end

function plot_infotransfer!(ax, M::Matrix, S::Matrix, vars::Vector, opts::NamedTuple;
    metric::String = "Normalized information transfer")
    heatmap!(ax, M, colormap = :balance, colorrange = (-100, 100))
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