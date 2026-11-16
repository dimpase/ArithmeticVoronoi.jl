# Changelog

## 0.1.1

- Construct certified real-quadratic SL₂ complexes through `hap_voronoi_complex(D; group=:SL)` and HAP's Julia provider. The determinant-height construction computes fresh cells and exact link-collapse certificates, with validation for discriminants 5 and 12. No Magma or C++ backend is required.
- Support rational GL/SL/PSL requests and full scalar-centre PSL quotients in the HAP interface, with direct Julia calls that need no manual provider registration.
- Keep the original standalone real-SL entry point as a compatibility wrapper. Detailed construction output and certificates remain available through `ArithmeticVoronoi.RealQuadraticSL`.

Real-quadratic higher-rank SL, real-quadratic PSL and nonfree-lattice complexes remain unsupported. A field whose chosen height cannot be certified produces an error; no stored result is substituted.

## 0.1.0

Initial registered release of the native Julia engine, with perfect forms, presentations, matrix-to-word computations and HAP integration.
