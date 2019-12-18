FMMPotential
============

FMMPotential is a wrapper around [exafmm][1], for calculating potentials and force vectors
in an electrostatic field.

The library bundles the content of the "3d" subfolder of a slightly patched exafmm version,
that was forked on 23. Aug. 2017.

The main changes in the exafmm code are:

- Using a static `std::vector` variables instead of variable length arrays.
- The tree is built iterative instead of recursive to avoid stack overflows when many points are coplanar.
- An index property was added to exafmm::Body type.
- Files not used by the FMMPotential library were removed.

[1]: https://github.com/exafmm/exafmm/

Compilation with OpenMP
-----------------------

The library must be compiled with a compiler that supports OpenMP 3.0, like clang or gcc.

### Using Visual Studio

To build the library with Microsoft Visual Studio, use the `LLVM-vs2014` toolchain, i.e., by using the
`-T"LLVM-vs2014"` option on the commandline or using "LLVM-vs2014" as the toolchain
option in cmake-gui with the MSVC project generator.  
See this [StackOverflow Answer](https://stackoverflow.com/a/38174328) for more details.

### Compiler flags for OpenMP

You need to make sure that OpenMP is found and the right libraries are linked in.
CMake does not find the default installation of clang in FindOpenMP on Windows automatically.

- When OpenMP is not found, add `-fopenmp` to compiler flags and the path to `libomp.lib`
(e.g. `C:\Program Files\LLVM\lib\libomp.lib`) to `CMAKE_CXX_STANDARD_LIBRARIES`.
- When using clang, use `-Xclang -fopenmp` in `OpenMP_CXX_FLAGS` and `OpenMP_C_FLAGS`.

### Runtime environment variable for multiple OpenMP libraries

When linking two different OpenMP libaries, like `libiomp5md.dll` / `libomp.lib` from clang
and `libomp` from MSVC, the environment variable `KMP_DUPLICATE_LIB_OK` must be set
to `TRUE` at runtime.
