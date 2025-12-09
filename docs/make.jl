using Documenter
using NCTSSOS

makedocs(
sitename = "NCTSSOS",
pages = ["Home" => "index.md",
         "Noncommutative Polynomial Optimization" => "ncpop.md",
         "Trace Polynomial Optimization" => "trace.md",
         "State Polynomial Optimization" => "state.md",
         "Sum-Of-Hermitian-Squares Optimization" => "sohs.md",
         "Examples" => ["Bell inequalities" => "bell.md",
         "Inclusion" => "inclusion.md",
         ],
         ],
modules = [NCTSSOS],
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/wangjie212/NCTSSOS.git"
)
