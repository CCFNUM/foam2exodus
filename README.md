# foam2exodus

A converter utility that translates OpenFOAM mesh and field data into the Exodus II database format.

## Overview

foam2exodus is a tool designed to convert OpenFOAM computational fluid dynamics (CFD) cases into the Exodus II format, a widely-used database format for finite element analysis. The Exodus II format enables visualization and post-processing in various scientific visualization tools such as ParaView, VisIt, and others.

### Current Status (Under Development)

**Currently Implemented:**
- Serial OpenFOAM mesh conversion (ASCII and binary formats)
- Conversion of mesh geometry (points, cells, faces)
- Boundary patch export as Exodus sidesets
- Support for hexahedral, tetrahedral, and pyramid cell types
- Cell zones preservation
- Conservative conversion of scalar, vector, spherical-tensor,
  symmetric-tensor and tensor volume fields to Exodus nodal variables
- Selection of a named OpenFOAM time or the latest time

**Planned Features:**
- Multi-time Exodus databases
- Parallel OpenFOAM case support
- Option to generate single or decomposed Exodus database files according to user preference

## Prerequisites

### Required Software

1. **CMake** (version 3.10 or higher)
   ```bash
   cmake --version  # Check your version
   ```

2. **C++ Compiler** with C++17 support
   - GCC 7.0+ or
   - Clang 5.0+ or
   - Any C++17 compliant compiler

3. **NetCDF Library** (with C interface)
   - Development headers and libraries are required

   Installation on Ubuntu/Debian:
   ```bash
   sudo apt-get install libnetcdf-dev pkg-config
   ```

   Installation on Fedora/RHEL/CentOS:
   ```bash
   sudo yum install netcdf-devel pkgconfig
   ```

   Installation on macOS (using Homebrew):
   ```bash
   brew install netcdf pkg-config
   ```

4. **pkg-config**
   - Used by CMake to locate NetCDF installation
   - Usually comes pre-installed on Linux systems
   - On macOS: `brew install pkg-config`

## Building from Source

### Step 1: Clone or Navigate to the Repository

```bash
cd /path/to/foam2exodus
```

### Step 2: Create Build Directory

```bash
mkdir build
cd build
```

### Step 3: Configure with CMake

```bash
cmake ..
```

If NetCDF is installed in a non-standard location, you may need to help CMake find it:
```bash
export PKG_CONFIG_PATH=/path/to/netcdf/lib/pkgconfig:$PKG_CONFIG_PATH
cmake ..
```

### Step 4: Compile

```bash
make
```

For faster compilation on multi-core systems:
```bash
make -j$(nproc)  # Linux
make -j$(sysctl -n hw.ncpu)  # macOS
```

### Step 5: Install (Optional)

```bash
sudo make install
```

This will install the `foam2exodus` executable to `/usr/local/bin` by default.

## Usage

### Basic Syntax

```bash
foam2exodus [--fields all|U,p,...] \
    [--field-names velocity,pressure,...] [--time latest|TIME] \
    <OpenFOAM_case_dir> <output.exo>
```

### Advanced

Convert multiple meshes and merge into single exodus file
```bash
foam2exodus --multiple <OpenFOAM_case_dir_1> <OpenFOAM_case_dir_2> <output.exo>
```

Convert selected fields from time `100`:
```bash
foam2exodus --time 100 --fields U,p,k /path/to/case solution.e
```

Convert selected fields and choose their Exodus base names:
```bash
foam2exodus --time 100 \
    --fields U,p,k \
    --field-names velocity,pressure,tke \
    /path/to/case solution.e
```

Convert all supported volume fields from the latest numeric time directory:
```bash
foam2exodus --fields all /path/to/case solution.e
```

### Arguments

- `<OpenFOAM_case_dir>`: Path to the OpenFOAM case directory containing the `constant/polyMesh` folder
- `<output.exo>`: Desired path and filename for the output Exodus II database file
- `--fields all`: Select every supported volume field in the chosen time
- `--fields U,p,k`: Select a comma-separated list (a quoted space-separated
  list such as `--fields '(U p k)'` is also accepted)
- `--field-names velocity,pressure,tke`: Rename the selected fields in the same
  order. There must be exactly one output base name per explicitly selected
  field. This option is intentionally rejected with `--fields all`.
- `--time latest|TIME`: Select the latest numeric time (the default when
  fields are requested) or an exact time-directory name

### Examples

Convert a serial OpenFOAM case:
```bash
./foam2exodus /path/to/openfoam/case myMesh.exo
```

Convert and specify full output path:
```bash
./foam2exodus ~/simulations/cavity ~/results/cavity_mesh.exo
```

### Output

The converter will display:
- Mesh statistics (number of points, cells, boundary patches)
- Progress information
- Success/error messages

Boundary patches from OpenFOAM are automatically exported as sidesets in the Exodus file, preserving their names and topology.

### Conservative CVFEM field projection

OpenFOAM volume fields are cell averages, whereas an Exodus nodal value used by
a CVFEM solver should represent the average over the median-dual control volume
surrounding that node. The converter therefore does not copy or inverse-distance
average cell-centre values.

For each OpenFOAM cell it:

1. Builds a weighted least-squares linear reconstruction from adjacent cell
   averages and values stored at physical boundary-face centroids.
2. Limits the reconstructed change on the stencil while retaining the original
   cell average.
3. Geometrically partitions each output element into median-dual
   sub-control-volume pieces and integrates the reconstruction over each piece.
4. Divides the accumulated integral at a node by that node's complete dual
   volume.

Consequently, a physical boundary-face value influences the adjacent boundary
nodes but never overwrites them. A boundary node is the volume average over its
whole dual control volume, including contributions from its adjacent cell
interiors. The projection preserves the global volume integral to roundoff and
prints the relative conservation residual for every Exodus variable.

`cyclic`, `cyclicAMI`, `processor`, and `empty` patches are not treated as
physical boundary samples. A physical patch entry without an explicit `value`
uses the owner-cell value, which is the conservative zero-normal-change
fallback.

Exodus stores scalar nodal variables, so components use these names:

- scalar: `p`
- vector: `U_x`, `U_y`, `U_z`
- symmetric tensor: `_xx`, `_xy`, `_xz`, `_yy`, `_yz`, `_zz`
- tensor: `_xx`, `_xy`, `_xz`, `_yx`, `_yy`, `_yz`, `_zx`, `_zy`, `_zz`

An output base name is applied before component suffixes. For example,
`--fields U,p --field-names velocity,pressure` produces `velocity_x`,
`velocity_y`, `velocity_z`, and `pressure`.

The selected OpenFOAM time is written to `time_whole`. Only volume fields are
projected; face fields such as `surfaceScalarField` are intentionally excluded
from `--fields all` because they do not represent element averages.

## OpenFOAM Mesh Format Support

### Supported Formats
- ASCII format (human-readable text files)
- Binary format (OpenFOAM binary mesh files)

### Supported Cell Types
- Hexahedra (hex)
- Tetrahedra (tet)
- Pyramids (pyramid)
- Prisms/wedges (wedge)
- Arbitrary polyhedra (e.g. cfMesh hex cells with hanging nodes): split into
  conformal tetrahedra/pyramids via a per-cell centroid and shared per-face
  centroids, so the result is watertight and every element has positive volume.
  A standard cell is decomposed the same way only if it fails validation.

### Cell Recognition and Validation

A cell is written as a standard Exodus element when it matches that element's
topology, which is tested the same way OpenFOAM's own `hexMatcher` and friends
do it: face count, face sizes, vertex count and how many faces each vertex is
shared by. Face counts alone are not enough.

The canonical node ordering comes from OpenFOAM's owner/neighbour convention
(a face normal points out of its owner cell), never from the cell geometry, so
slivers and strongly warped cells are ordered as reliably as regular ones. That
convention yields exactly the ordering OpenFOAM's `cellShape` models use, which
is also the Exodus node ordering.

Before anything is written, every element is checked for supported topology,
node count, repeated node IDs, coincident vertices, collapsed edges, corner
Jacobian signs and signed volume. HEX8 volumes use the trilinear (2x2x2 Gauss)
formula, the same one ParaView reports. A cell that fails is named in the
message by its OpenFOAM cell ID; if no valid element can be produced at all the
conversion fails rather than writing an inverted or collapsed element.

### Mesh Components
- Points (vertices)
- Faces
- Cell connectivity
- Boundary patches
- Cell zones
- Owner/neighbour cell relationships

## Project Structure

```
foam2exodus/
├── CMakeLists.txt           # CMake build configuration
├── README.md                # This file
├── foam2exodus/
│   ├── main.cpp            # Main entry point
│   ├── src/
│   │   ├── OpenFOAMMeshReader.h    # OpenFOAM mesh reader interface
│   │   ├── OpenFOAMMeshReader.cpp  # OpenFOAM mesh parser implementation
│   │   ├── ExodusWriter.h          # Exodus II writer interface
│   │   └── ExodusWriter.cpp        # Exodus II database writer
│   ├── examples/           # Example test cases
│   └── tests/              # Regression tests and their polyMesh fixtures
└── build/                  # Build directory (created by user)
```

## Tests

```bash
python3 tests/run_tests.py            # uses build/foam2exodus
python3 tests/run_tests.py --exe /path/to/foam2exodus
```

The fixtures in `tests/meshes/` cover an all-hexahedral mesh (a 280-cell cut-out
of the WILO impeller mesh around its worst high-aspect-ratio cells), a genuinely
polyhedral mesh (the dual of a block mesh), and the mixed hex/tet `3DCube`
example. The tests check element counts and types, connectivity, orientation,
volumes, scaled Jacobians, side-set coverage and source-cell traceability.
Requires `numpy` and either `netCDF4` or `scipy`.

## Troubleshooting

### CMake Cannot Find NetCDF

If you encounter:
```
Could not find module FindNetCDF.cmake
```

Ensure:
1. NetCDF development packages are installed
2. pkg-config is installed and functional
3. PKG_CONFIG_PATH includes the NetCDF pkgconfig directory

Test NetCDF installation:
```bash
pkg-config --modversion netcdf
```

### Compilation Errors

If you see C++17-related errors, ensure your compiler supports C++17:
```bash
g++ --version  # Should be 7.0 or higher
clang++ --version  # Should be 5.0 or higher
```

### Runtime Errors

**"Cannot find polyMesh directory"**
- Ensure the OpenFOAM case directory contains `constant/polyMesh/`
- Verify the case path is correct

**"Error reading mesh files"**
- Check that mesh files (points, faces, owner, neighbour, boundary) exist
- Verify file permissions are readable

## Contributing

This project is under active development. Contributions, bug reports, and feature requests are welcome.

## Licence

BSD-3-Clause.  See file headers for details.

Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and Arts.

## Contact

* **Project Maintainer:** Lucian Hanimann ([lucian.hanimann@hslu.ch](mailto:lucian.hanimann@hslu.ch))
