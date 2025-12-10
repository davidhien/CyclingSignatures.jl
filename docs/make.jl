using CyclingSignatures
using Documenter

DocMeta.setdocmeta!(CyclingSignatures, :DocTestSetup, :(using CyclingSignatures); recursive=true)

makedocs(;
    modules=[CyclingSignatures],
    authors="davidhien <david.hien@outlook.de> and contributors",
    sitename="CyclingSignatures.jl",
    format=Documenter.HTML(;
        canonical="https://davidhien.github.io/CyclingSignatures.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "General Usage" => "general-usage.md",
        "API" => "api.md",
        "Internal" => [
            "CyclingSignatures" => "cycling-signatures.md",
            "ATTools" => "at-tools.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/davidhien/CyclingSignatures.jl",
    devbranch="main",
)
