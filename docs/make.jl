using Documenter
using NCTSSOS

makedocs(;
    sitename="NCTSSOS",
    pages=[
        "Home" => "index.md",
        "API" => "api.md"
    ],
    modules=[NCTSSOS],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
)

deploydocs(; repo="github.com/wangjie212/NCTSSOS.git")
