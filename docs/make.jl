using RareLotkaVolterra
using Documenter

DocMeta.setdocmeta!(
    RareLotkaVolterra,
    :DocTestSetup,
    :(using RareLotkaVolterra);
    recursive = true,
)

# GitLab repository where the package is hosted.
repo = "https://gitlab.com/ismael-lajaaiti/RareLotkaVolterra.jl/"

makedocs(;
    modules = [RareLotkaVolterra],
    authors = "Ismaël Lajaaiti <ismael.lajaaiti@gmail.com> and contributors",
    repo = repo * "blob/{commit}{path}#{line}",
    sitename = "RareLotkaVolterra.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)
