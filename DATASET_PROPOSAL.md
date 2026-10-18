# Proposal for reusable arithmetic-group data

This local dataset provides exact presentations and constructive matrix words
for GL₂, SL₂ and PSL₂ over real-quadratic maximal orders. The initial fields
are the four Kirchheimer–Wolfart discriminants 5, 8, 12 and 13, followed by
the ordinary-class-number-two test D40. Each advertised completion status is
backed by the independent Voronoi validation and coverage data; unfinished
stages remain separately labelled.

A possible LMFDB contribution would link a compact arithmetic-group record
to its number-field record and provide downloadable presentations, exact
matrix realizations, abelianizations, generator maps and constructive-word
examples. Stable arithmetic-object identifiers would allow several
presentations and computation versions for one group. No new LMFDB web
interface is proposed at this stage.

The large exact Voronoi/attaching-path certificates could live in a versioned
citable archive, with its manifest hash, software/environment and independent
verification commands linked from the compact record. The current archive
is a local candidate without a DOI. The schema distinguishes full GL, the
intermediate fundamental-scalar quotient, SL and PSL explicitly. It records
ordinary class groups, not an unlabelled mixture of ordinary and narrow data.

The contribution would add group-presentation data alongside existing
[number-field records](https://www.lmfdb.org/NumberField/) and related
[Hilbert modular form data](https://www.lmfdb.org/ModularForm/GL2/TotallyReal/).
Field labels and generator correspondences remain null until verified.
Literature comparisons for D5/D8/D12/D13 retain exact substitutions and
central signs, and identify the supplied D8 transcription discrepancy.

Questions for a later editorial discussion are the desired group-table scope,
field linkage conventions, presentation versioning, certificate-archive
hosting and the level of verification editors would want to rerun. No editor
has been contacted and no dataset, issue or pull request has been submitted
as part of this task.
