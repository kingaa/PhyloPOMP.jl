module PhyloPOMPMakieExt

using PhyloPOMP
import PhyloPOMP: genealogyplot, treelayout
using Makie

genealogyplot(
    g::Genealogy;
    ladderize::Bool = true,
    points::Bool = true,
    palette = Makie.wong_colors(),
    title::AbstractString = "",
    figure = (;),
    axis = (;),
) = begin
    L = treelayout(g; ladderize)
    fig = Figure(; figure...)
    ax = Axis(
        fig[1,1];
        xlabel = "time", title = title,
        yticksvisible = false, yticklabelsvisible = false, ygridvisible = false,
        limits = (g.t0, g.time, nothing, nothing),
        axis...,
    )
    segments = Point2f[]
    for (x0, x1, y) ∈ L.branches
        push!(segments, Point2f(x0, y), Point2f(x1, y))
    end
    for (x, y0, y1) ∈ L.connectors
        push!(segments, Point2f(x, y0), Point2f(x, y1))
    end
    linesegments!(ax, segments; color = :black)
    if points
        s = Int.(samples(g))
        demes = unique(g[i].deme for i ∈ s if !ismissing(g[i].deme))
        if isempty(demes)
            scatter!(ax, [g[i].slate for i ∈ s], [L.y[i] for i ∈ s]; color = palette[1])
        else
            for d ∈ sort(demes; by = Int)
                j = [i for i ∈ s if isequal(g[i].deme, d)]
                scatter!(
                    ax, [g[i].slate for i ∈ j], [L.y[i] for i ∈ j];
                    color = palette[mod1(Int(d), length(palette))], label = string(d),
                )
            end
            axislegend(ax; position = :lt)
        end
    end
    fig
end

end # module
