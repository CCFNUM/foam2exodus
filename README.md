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

**Planned Features:**
- Field data conversion from all time directories
- Parallel OpenFOAM case support
- Option to generate single or decomposed Exodus database files according to user preference
- Time-series data export

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
foam2exodus <OpenFOAM_case_dir> <output.exo>
```

### Advanced

Convert multiple meshes and merge into single exodus file
```bash
foam2exodus --multiple <OpenFOAM_case_dir_1> <OpenFOAM_case_dir_2> <output.exo>
```

### Arguments

- `<OpenFOAM_case_dir>`: Path to the OpenFOAM case directory containing the `constant/polyMesh` folder
- `<output.exo>`: Desired path and filename for the output Exodus II database file

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

## OpenFOAM Mesh Format Support

### Supported Formats
- ASCII format (human-readable text files)
- Binary format (OpenFOAM binary mesh files)

### Supported Cell Types
- Hexahedra (hex)
- Tetrahedra (tet)
- Pyramids (pyramid)

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
│   └── examples/           # Example test cases
└── build/                  # Build directory (created by user)
```

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