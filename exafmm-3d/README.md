This code was extracted from the "3d" subfolder of the [exafmm](https://github.com/exafmm/exafmm)
project on August 23, 2017 and contains a few changes by the authors of the FMMPotential library.

The main changes are:

- Using a static `std::vector` variables instead of variable length arrays.
- The tree is built iterative instead of recursive to avoid stack overflows when many points are coplanar.
- An index property was added to exafmm::Body type.
- Files not used by the FMMPotential library were removed.
