using SpeciesReactivity
using Documenter

DocMeta.setdocmeta!(
    SpeciesReactivity,
    :DocTestSetup,
    :(using SpeciesReactivity);
    recursive = true,
)

# GitLab repository where the package is hosted.
repo = "https://github.com/ismael-lajaaiti/species-reactivity"

makedocs(;
    modules = [SpeciesReactivity],
    authors = "Ismaël Lajaaiti <ismael.lajaaiti@gmail.com> and contributors",
    repo = repo * "blob/{commit}{path}#{line}",
    sitename = "SpeciesReactivity.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)
