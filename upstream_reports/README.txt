Upstream reports submitted on 2026-10-17 as requested by the user
================================================================

Three separate issues were filed from the dimpase account:

  #13 PERF_SerialPerfectComputation ignores ComputeComplex=F
  https://github.com/MathieuDutSik/polyhedral_common/issues/13

  #14 Propagate failed boundary-product validation instead of returning the complex
  https://github.com/MathieuDutSik/polyhedral_common/issues/14

  #15 PERFCOMP_general_query: handle paths with spaces and nonzero child exits
  https://github.com/MathieuDutSik/polyhedral_common/issues/15

The ignored option and GAP process failures were reproduced. The boundary
report is explicitly a source-level finding; no actual nonzero boundary
product is claimed. These are issues, not code-fix PRs. No full internal audit
was posted and no upstream source or branches were changed for these reports.

The numbered .md files are the exact issue bodies. requests/ contains the
submitted JSON payloads; posted/ contains readbacks from GitHub. The submission
manifest verifies matching text, title, author and open status.

The commands below use the original build layout and record the submitted
reproductions. For current PerfectForms build instructions, see USER_GUIDE.md.
To reproduce these probes, supply the corresponding upstream checkout and
driver paths as arguments:

These optional audit probes also require a standalone GAP executable on PATH
(or --gap passed to reproduce.py). This is separate from the GAP dependency
installed by the Julia project for ordinary computations.

  cmake --build build --target PERF_SerialPerfectComputation --parallel 1
  python3 upstream_reports/reproduce.py \
    --source build/_deps/polyhedral_common-src \
    --binary build/bin/PERF_SerialPerfectComputation

The standalone GAP script embedded in #15 was also executed successfully:

  python3 upstream_reports/reproduce_gap.py build/_deps/polyhedral_common-src

reproduction_summary.json identifies the accepted fresh logs/input directories.
Earlier development failures remain as diagnostics and are not the submitted
reproduction results. The original mathematical implementation is unchanged.

Generated reproduction directories, including extracted upstream query
functions, are now ignored local diagnostics rather than published files.
The summary, submitted reports, receipts and our reproduction scripts remain
versioned. To reproduce a report, supply the separately fetched pinned
upstream checkout; the scripts read its source at execution time.
