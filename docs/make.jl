using DynamicPanel
using Documenter

DocMeta.setdocmeta!(DynamicPanel, :DocTestSetup, :(using DynamicPanel); recursive=true)

makedocs(;
    modules=[DynamicPanel],
    authors="Raymond Chiasson",
    repo="https://github.com/RC416/DynamicPanel.jl/blob/{commit}{path}#{line}",
    sitename="DynamicPanel.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RC416.github.io/DynamicPanel.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/RC416/DynamicPanel.jl",
    devbranch="master",
)
