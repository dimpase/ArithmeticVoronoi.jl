`PERF_SerialPerfectComputation` accepts `ComputeComplex = F`, but still constructs the full complex. Could this option either control construction, or be removed/documented as having no effect on this driver?

Tested at `cfcd79dbddf67a494a1c71b776174b190da5209a`, currently `master`, using the GMP driver. This came up while using the real-quadratic perfect-form enumeration from OSCAR: we need perfect forms and labelled neighbours, without all lower-dimensional faces.

### Reproduction

In a fresh working directory, save as `d5.nml`:

```text
&TSPACE
 TypeTspace = "RealQuad"
 ClassicDim = 2
 RealImagDim = 2
 RealImagSum = 1
 RealImagProd = -1
 SuperMatMethod = "NotNeeded"
 ListComm = "Use_realimag"
 PtGroupMethod = "Trivial"
 FileListSubspaces = "unset"
/
&DATA
 arithmetic = "gmp"
 OnlyWellRounded = T
 ComputeComplex = F
 ComputeBoundary = F
 ComputePreImage = F
 CacheFile = "none"
/
&QUERIES
 FileListNumberOrbit = "orbit_counts.g"
/
```

Run the built executable, adjusting its absolute path:

```sh
timeout 45s prlimit --as=10737418240 -- /absolute/path/to/PERF_SerialPerfectComputation d5.nml
cat orbit_counts.g
```

Observed, with exit status 0:

```gap
return [2,4,5,4,1,0];
```

Repeat in a separate fresh directory with only `ComputeComplex = T` changed. The output is identical. Thus `F` does not restrict construction to the perfect forms at level zero; all lower levels are present. We also observed this with the classical rank-two example (`[1,1,0]`). This is not a timeout report: both D5 runs complete in about 0.115 s locally.

### Source

- The option is registered, defaulting to false, in [perfect_tspace.h](https://github.com/MathieuDutSik/polyhedral_common/blob/cfcd79dbddf67a494a1c71b776174b190da5209a/src_perfect/perfect_tspace.h#L411-L430).
- [process_A](https://github.com/MathieuDutSik/polyhedral_common/blob/cfcd79dbddf67a494a1c71b776174b190da5209a/src_perfect/PERF_SerialPerfectComputation.cpp#L37-L50) reads the other complex options but not `ComputeComplex`, and always calls `get_full_complex_enumeration`.

The intended semantics need clarification because some queries inherently need the full complex. A perfect-form-only entry point would also solve our integration need. If the option is obsolete, rejecting it or documenting that this executable always builds all levels would avoid suggesting that it limits the work.

For context, the underlying D5 perfect-form data successfully passed independent exact OSCAR/Hecke checks: two GL(2,O) orbits, stabilizer orders 20 and 12, and eight labelled facet orbits. This report concerns option handling, not an observed arithmetic failure.

Build: GCC 16.2.0, C++20, GMP 6.3.0, Eigen 3.4.1, Boost 1.90, `-O1 -DNDEBUG -DSANITY_CHECK`; serial execution. The nauty fork was pinned to `cb9dc822f9d1aa3a02b1dcb7a18dd9609d838b1e`. The driver source and mathematical headers were unmodified.
