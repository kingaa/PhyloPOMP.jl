using Documenter
using PhyloPOMP

makedocs(
    sitename = "PhyloPOMP.jl",
    modules  = [PhyloPOMP],
    repo = Remotes.GitHub("kingaa","PhyloPOMP.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true", # easier local build
        size_threshold = 600 * 2^10,
        size_threshold_warn = 500 * 2^10,
        canonical="https://kingaa.github.io/PhyloPOMP.jl/stable/",
        edit_link=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => ["seir.md"],
    ]
)

deploydocs(
    ;
    repo="github.com/kingaa/PhyloPOMP.jl"
)
