# Release notes

## 0.1.0

First release proposed for Julia's General registry. Computes perfect forms,
GL(n,O) presentations and matrix words for free lattices over maximal
quadratic orders, as well as GL(n,Z); the rank-two pipeline also computes
SL and PSL presentations. Includes a native provider for HAP cell complexes.

Installation uses registered Julia dependencies. Computation defaults to the
native Julia/Hecke/OSCAR backend, with output under the caller's working
directory. The C++ comparison adapter uses a separate, explicitly loaded
environment. Recorded presentations remain in the optional
ArithmeticVoronoiData package and are never loaded implicitly by the engine.
