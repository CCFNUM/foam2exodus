cmake_minimum_required(VERSION 3.10)
project(OpenFOAMToExodus)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED netcdf)

include_directories(${NETCDF_INCLUDE_DIRS})


# Add src directory to include path
include_directories(${CMAKE_SOURCE_DIR}/src)

# Add executable with source files from src/
add_executable(foam2Exodus
    main.cpp
    src/OpenFOAMMeshReader.cpp
    src/ExodusWriter.cpp
)

target_link_libraries(foam2Exodus ${NETCDF_LIBRARIES})

install(TARGETS foam2Exodus DESTINATION bin)
