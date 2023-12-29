using Documenter
using JuliaFormatter
using Random
using SpeciesReactivity
using Test

seed = rand(1:10^6)
Random.seed!(seed)
@info "Seed set to $seed."

DocMeta.setdocmeta!(
    SpeciesReactivity,
    :DocTestSetup,
    :(using SpeciesReactivity);
    recursive = true,
)
doctest(SpeciesReactivity)
println("------------------------------------------")

test_files = [
# Add test files here.
]

# Text formatting shorthand.
highlight = "\033[7m"
bold = "\033[1m"
green = "\033[32m"
reset = "\033[0m"

no_break = true
for test in test_files
    println("$(highlight)$(test)$(reset)")
    global no_break = false
    include(test) # If a test fails, the loop is broken.
    global no_break = true
    println("$(bold)$(green)PASSED$(reset)")
    println("------------------------------------------")
end

if no_break
    @info "Checking source code formatting.."
    exclude = [
        "CONTRIBUTING.md", # Not formatted according to JuliaFormatter.
    ]
    for (folder, _, files) in walkdir("..")
        for file in files
            path = joinpath(folder, file)
            display_path = joinpath(splitpath(path)[2:end]...)
            if display_path in exclude
                continue
            end
            if !any(endswith(file, ext) for ext in [".jl", ".md", ".jmd", ".qmd"])
                continue
            end
            println(display_path)
            if !format(path; overwrite = false, format_markdown = true)
                config_path =
                    joinpath(basename(dirname(abspath(".."))), ".JuliaFormatter.toml")
                dev_path = escape_string(abspath(path))
                @warn "Source code in $file is not formatted according \
                to the project style defined in $config_path. \
                Consider formatting it using your editor's autoformatter or with:\n\
                    julia> using JuliaFormatter;\n\
                    julia> format(\"$dev_path\", format_markdown=true)\n"
            end
        end
    end
end
