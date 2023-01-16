using Documenter
using NCTSSOS

makedocs(
sitename = "NCTSSOS",
pages = ["Home" => "index.md",
         "Noncommutative Polynomial Optimization" => "ncpop.md"],
modules = [NCTSSOS],
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/wangjie212/NCTSSOS.git"
)
