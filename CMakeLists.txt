cmake_minimum_required(VERSION 3.10)
project(OpenFOAMToExodus VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Portability settings
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(PkgConfig REQUIRED)

# NetCDF configuration
if(DEFINED NETCDF_DIR)
    # Use custom NetCDF path
    message(STATUS "Using custom NetCDF path: ${NETCDF_DIR}")
    set(NETCDF_INCLUDE_DIRS "${NETCDF_DIR}/include")
    set(NETCDF_LIBRARY_DIRS "${NETCDF_DIR}/lib")
    set(NETCDF_BIN_DIR "${NETCDF_DIR}/bin")

    # Find the netcdf library
    find_library(NETCDF_LIBRARY
        NAMES netcdf
        PATHS ${NETCDF_LIBRARY_DIRS}
        NO_DEFAULT_PATH
    )

    if(NETCDF_LIBRARY)
        set(NETCDF_LIBRARIES ${NETCDF_LIBRARY})
        set(NETCDF_FOUND TRUE)
        link_directories(${NETCDF_LIBRARY_DIRS})
        message(STATUS "Found NetCDF library: ${NETCDF_LIBRARY}")
        message(STATUS "NetCDF include dir: ${NETCDF_INCLUDE_DIRS}")
        message(STATUS "NetCDF library dir: ${NETCDF_LIBRARY_DIRS}")
        message(STATUS "NetCDF bin dir: ${NETCDF_BIN_DIR}")
    else()
        message(FATAL_ERROR "NetCDF library not found in ${NETCDF_DIR}/lib")
    endif()
else()
    # Use pkg-config to find NetCDF (standard behavior)
    pkg_check_modules(NETCDF REQUIRED netcdf)
endif()

include_directories(${NETCDF_INCLUDE_DIRS})


# Add src directory to include path
include_directories(${CMAKE_SOURCE_DIR}/src)

# CLI executable
add_executable(foam2Exodus
    main.cpp
    src/OpenFOAMMeshReader.cpp
    src/ExodusWriter.cpp
    src/MergedMeshReader.cpp
)

target_link_libraries(foam2Exodus ${NETCDF_LIBRARIES})

install(TARGETS foam2Exodus DESTINATION bin)

# GUI executable (optional, only if Qt5 and GUI sources are found)
if(DEFINED QT5_DIR)
    # Use custom Qt5 path
    message(STATUS "Using custom Qt5 path: ${QT5_DIR}")
    set(CMAKE_PREFIX_PATH "${QT5_DIR}")
    find_package(Qt5 COMPONENTS Widgets QUIET)
else()
    # Use standard Qt5 detection
    find_package(Qt5 COMPONENTS Widgets QUIET)
endif()

if(Qt5_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/src/main_gui.cpp")
    message(STATUS "Qt5 and GUI sources found, building GUI application")

    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)
    set(CMAKE_AUTOUIC ON)

    add_executable(foam2ExodusGUI
        src/main_gui.cpp
        src/MainWindow.cpp
        src/OpenFOAMMeshReader.cpp
        src/ExodusWriter.cpp
        src/MergedMeshReader.cpp)

    target_link_libraries(foam2ExodusGUI
        ${NETCDF_LIBRARIES}
        Qt5::Widgets)

    install(TARGETS foam2ExodusGUI DESTINATION bin)
else()
    if(NOT Qt5_FOUND)
        message(STATUS "Qt5 not found, skipping GUI application build")
        message(STATUS "Install Qt5 development packages to build the GUI:")
        message(STATUS "  Ubuntu/Debian: sudo apt install qtbase5-dev")
        message(STATUS "  Fedora/RHEL: sudo dnf install qt5-qtbase-devel")
        message(STATUS "  Arch: sudo pacman -S qt5-base")
    else()
        message(STATUS "GUI source files not found, skipping GUI application build")
    endif()
endif()
