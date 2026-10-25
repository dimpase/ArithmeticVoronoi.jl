# Optional integration suite; requires a built and configured PerfectForms
# exporter. See USER_GUIDE.md for explicit historical-fixture restoration.
include(joinpath(@__DIR__,"../integrations/perfect_forms/adapter.jl"))
using .PerfectFormsAdapter
include("checkpoint_opt_in.jl")
include("contracts.jl")
include("two_cells.jl")
include("fields_words.jl")
FIELD_WORD_RESULTS=Dict(D=>run_field_pipeline(D;backend=:polyhedral_common) for D in (5,8,12,13,40))
run_fields_words_tests(FIELD_WORD_RESULTS)
