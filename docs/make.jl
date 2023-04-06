using RareLotkaVolterra
using Documenter

DocMeta.setdocmeta!(RareLotkaVolterra, :DocTestSetup, :(using RareLotkaVolterra); recursive=true)

makedocs(;
    modules=[RareLotkaVolterra],
    authors="Ismaël Lajaaiti <ismael.lajaaiti@gmail.com> and contributors",
    repo="https://gitlab.com/Ismaël Lajaaiti/RareLotkaVolterra.jl/blob/{commit}{path}#{line}",
    sitename="RareLotkaVolterra.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
