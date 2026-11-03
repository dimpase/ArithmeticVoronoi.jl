using ArithmeticVoronoi, Documenter

DocMeta.setdocmeta!(ArithmeticVoronoi, :DocTestSetup,
    :(using ArithmeticVoronoi, Oscar); recursive=true)

makedocs(
    modules=[ArithmeticVoronoi],
    sitename="ArithmeticVoronoi.jl",
    authors="Dima Pasechnik and contributors",
    checkdocs=:exports,
    format=Documenter.HTML(
        prettyurls=get(ENV,"CI","false")=="true",
        canonical="https://dimpase.github.io/ArithmeticVoronoi.jl/",
    ),
    pages=[
        "Home"=>"index.md",
        "User guide"=>"guide.md",
        "Computation requests"=>"requests.md",
        "Algorithms"=>[
            "Native enumeration"=>"native_backend.md",
            "Presentations"=>"theory/presentations.md",
            "Matrix words"=>"theory/words.md",
        ],
        "Examples and validation"=>[
            "Real quadratic fields"=>"native_survey.md",
            "Imaginary quadratic fields"=>"imaginary_quadratic.md",
            "Higher rank"=>"higher_rank.md",
            "HAP and Magma"=>"hap_and_vmh.md",
            "Kirchheimer–Wolfart comparison"=>"validation/kirchheimer_wolfart.md",
            "Historical C++ backend audit"=>"validation/backend_audit.md",
        ],
        "Data and reproducibility"=>[
            "Code and data separation"=>"code_and_data.md",
            "Archive schema"=>"data/schema.md",
            "Data library proposal"=>"data/proposal.md",
            "Repository transition"=>"repository_transition.md",
        ],
        "API reference"=>"api.md",
        "Development"=>"development.md",
    ],
)
