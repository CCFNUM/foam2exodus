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
pkg_check_modules(NETCDF REQUIRED netcdf)

include_directories(${NETCDF_INCLUDE_DIRS})


# Add src directory to include path
include_directories(${CMAKE_SOURCE_DIR}/src)

# CLI executable
add_executable(foam2Exodus
    main.cpp
    src/OpenFOAMMeshReader.cpp
    src/ExodusWriter.cpp
)

target_link_libraries(foam2Exodus ${NETCDF_LIBRARIES})

install(TARGETS foam2Exodus DESTINATION bin)

# GUI executable (optional, only if Qt5 is found)
find_package(Qt5 COMPONENTS Widgets QUIET)

if(Qt5_FOUND)
    message(STATUS "Qt5 found, building GUI application")

    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)
    set(CMAKE_AUTOUIC ON)

    add_executable(foam2ExodusGUI
        src/main_gui.cpp
        src/MainWindow.cpp
        src/OpenFOAMMeshReader.cpp
        src/ExodusWriter.cpp)

    target_link_libraries(foam2ExodusGUI
        ${NETCDF_LIBRARIES}
        Qt5::Widgets)

    install(TARGETS foam2ExodusGUI DESTINATION bin)
else()
    message(STATUS "Qt5 not found, skipping GUI application build")
    message(STATUS "Install Qt5 development packages to build the GUI:")
    message(STATUS "  Ubuntu/Debian: sudo apt install qtbase5-dev")
    message(STATUS "  Fedora/RHEL: sudo dnf install qt5-qtbase-devel")
    message(STATUS "  Arch: sudo pacman -S qt5-base")
endif()
