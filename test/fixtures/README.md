# Regression fixtures

These files are test inputs, never imported by the computation engine.

`cpp_D5_export.json` is an existing raw C++ export for discriminant 5,
produced with polyhedral_common commit
`88ba28d184c49b0f27b83d2dd3681e3adc928629`. The checkpoint opt-in test uses it
with a deliberately failing test executable to distinguish fresh execution
from explicitly requested replay. Replay still performs the full exact
arithmetic validation of the stored export. This control-flow test does not
replace the real-exporter integration tests in `test/contracts.jl`.

`two_cells_v1.tar.gz` and its manifest are historical geometric certificates
used by the explicit regression-fixture restoration script. Production
pipelines do not automatically select those archived certificates.
