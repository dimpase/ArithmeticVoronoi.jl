# Changelog

## 0.1.2

- Generalize certified real-quadratic SL complexes to every rank `n ≥ 2`, retaining the correct scalar centre in odd and even rank.
- Fix the stabilizer check in equivariant cancellation: for a boundary term g·f of t, compare Stab(t) with g Stab(f) g⁻¹. A regression test uses noncentral stabilizers.
- Verify H₁ = 0 for SL₃ over the integers of Q(√5) against an independent presentation and a Magma/HAP calculation.
- Construct GL boundaries from local face incidences and certify Morse links by exact graph collapses, avoiding exhaustive flag expansion.
- Compute the cell skeleton needed by a requested HAP resolution and report its dimension explicitly. The detailed construction interface still offers the full geometric level.
- Generalize the independent Magma/HAP comparison to the index-`2n` subgroup `SL_n × Z`; reference data remain outside the engine.
- Keep certified ball sections as cells, subdivide other sections compatibly, and reduce the required chain modules with exact orientation and boundary checks.
- Avoid HAP coset searches for proved trivial orientation characters; retain the general formula for all other cases.

## 0.1.1

- Construct certified real-quadratic SL₂ complexes through `hap_voronoi_complex(D; group=:SL)` and HAP's Julia provider. The determinant-height construction computes fresh cells and exact link-collapse certificates, with validation for discriminants 5 and 12. No Magma or C++ backend is required.
- Support rational GL/SL/PSL requests and full scalar-centre PSL quotients in the HAP interface, with direct Julia calls that need no manual provider registration.
- Keep the original standalone real-SL entry point as a compatibility wrapper. Detailed construction output and certificates remain available through `ArithmeticVoronoi.RealQuadraticSL`.

Real-quadratic higher-rank SL, real-quadratic PSL and nonfree-lattice complexes remain unsupported. A field whose chosen height cannot be certified produces an error; no stored result is substituted.

## 0.1.0

Initial registered release of the native Julia engine, with perfect forms, presentations, matrix-to-word computations and HAP integration.
