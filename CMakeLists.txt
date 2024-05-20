# Defines the CMake commands/policies
cmake_minimum_required(VERSION 3.15)


# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Setup the project
project(geouned)
# target_compile_features(my_executable PRIVATE cxx_std_11)

include(GNUInstallDirs)



# Include the FetchContent module
include(FetchContent)

# Declare the fmt library
# FetchContent_Declare(
#   fmt
#   GIT_REPOSITORY https://github.com/fmtlib/fmt.git
#   GIT_TAG 8.0.1
# )

# # Make fmt available to the rest of your project
# FetchContent_MakeAvailable(fmt)

FetchContent_Declare(
  xerces-c
  GIT_REPOSITORY https://github.com/apache/xerces-c.git
  GIT_TAG v3.2.3
)

FetchContent_GetProperties(xerces-c)
if(NOT xerces-c_POPULATED)
  FetchContent_Populate(xerces-c)
  add_subdirectory(${xerces-c_SOURCE_DIR} ${xerces-c_BINARY_DIR})
  list(APPEND CMAKE_PREFIX_PATH ${xerces-c_BINARY_DIR})
endif()

FetchContent_MakeAvailable(xerces-c)

FetchContent_Declare(
  libmed
  GIT_REPOSITORY https://github.com/med-file/med.git
  GIT_TAG v4.1.0
)

FetchContent_GetProperties(libmed)
if(NOT libmed_POPULATED)
  FetchContent_Populate(libmed)
  add_subdirectory(${libmed_SOURCE_DIR} ${libmed_BINARY_DIR})
  list(APPEND CMAKE_PREFIX_PATH ${libmed_BINARY_DIR})
endif()

FetchContent_MakeAvailable(libmed)

include(ExternalProject)

ExternalProject_Add(freecad-project
  GIT_REPOSITORY https://github.com/FreeCAD/FreeCAD.git
  GIT_TAG 0.21.2
  GIT_SUBMODULES ""
  SOURCE_DIR ${CMAKE_BINARY_DIR}/freecad-source
  BINARY_DIR ${CMAKE_BINARY_DIR}/freecad-build
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -S ${CMAKE_BINARY_DIR}/freecad-source
    -B ${CMAKE_BINARY_DIR}/freecad-build
    -DBUILD_DOC=OFF
    -DFREECAD_USE_EXTERNAL_SMESH=OFF
    -DBUILD_FLAT_MESH=OFF
    -DBUILD_FEM_NETGEN=OFF
    -DXercesC_INCLUDE_DIR=${xerces-c_SOURCE_DIR}/src
    -DXercesC_LIBRARY=${xerces-c_BINARY_DIR}/src/XercesC
    -DXercesC_VERSION=3.2.3
    -DFREECAD_USE_EXTERNAL_SMESH=OFF
    -DMEDFILE_INCLUDE_DIR=${libmed_SOURCE_DIR}/include
    -DMEDFILE_LIBRARY=${libmed_BINARY_DIR}/lib
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/freecad-build
  INSTALL_COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR}/freecad-build
)


# cmake_minimum_required(VERSION 3.14)
# project(FreeCAD)

# include(FetchContent)

# FetchContent_Declare(
#   freecad
#   GIT_REPOSITORY https://github.com/FreeCAD/FreeCAD.git
#   GIT_TAG 0.21.2
#   GIT_SUBMODULES ""  # Fetch all submodules
# )

# FetchContent_GetProperties(freecad)
# if(NOT freecad_POPULATED)
#   FetchContent_Populate(freecad)
#   add_subdirectory(${freecad_SOURCE_DIR} ${freecad_BINARY_DIR})
# endif()

# # Find Python interpreter, libraries, etc.
# find_package(Python3 COMPONENTS Interpreter Development)

# # Find Python's setuptools package
# find_package(Python3 COMPONENTS setuptools)

# # Add the FreeCAD Python module
# Python3_add_library(freecad MODULE ${freecad_SOURCE_DIR}/src/Mod/FreeCAD.py)

# # Install the FreeCAD Python module
# install(TARGETS freecad
#   LIBRARY DESTINATION ${Python3_SITELIB}
# )