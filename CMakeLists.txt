# prepair
#
# Copyright © 2009-2026,
# Ken Arroyo Ohori    k.ohori@tudelft.nl
# Hugo Ledoux         h.ledoux@tudelft.nl
# Martijn Meijers     b.m.meijers@tudelft.nl
# All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.20)

project(pprepair VERSION 1.1.0 LANGUAGES CXX)

# Use a sensible default for single-config generators. Multi-config generators
# (such as Visual Studio) choose the configuration at build time instead.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# Dependencies
find_package(CGAL 6 REQUIRED)
find_package(GDAL 3.0 REQUIRED)

add_executable(pprepair
  FaceInfo.cpp
  IOWorker.cpp
  PlanarPartition.cpp
  PolygonHandle.cpp
  pprepair.cpp
)

target_compile_features(pprepair PRIVATE cxx_std_17)
set_target_properties(pprepair PROPERTIES CXX_EXTENSIONS OFF)
target_compile_definitions(pprepair PRIVATE PPREPAIR_VERSION="${PROJECT_VERSION}")

if(MSVC)
  # CGAL's constrained triangulation templates exceed the default object
  # section limit on MSVC, and the sources contain UTF-8 comments.
  target_compile_options(pprepair PRIVATE /bigobj /utf-8)
endif()

if(WIN32)
  target_compile_definitions(pprepair PRIVATE
    NOMINMAX
    WIN32_LEAN_AND_MEAN
    _CRT_SECURE_NO_WARNINGS
    _USE_MATH_DEFINES
  )
endif()

target_link_libraries(pprepair PRIVATE CGAL::CGAL GDAL::GDAL)

include(GNUInstallDirs)
install(TARGETS pprepair RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

enable_testing()
add_test(NAME pprepair_help COMMAND pprepair --help)
