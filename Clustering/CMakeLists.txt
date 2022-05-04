cmake_minimum_required (VERSION 3.5)
project(pyClustering)
find_package(pybind11 REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(/Users/lg9783/miniconda3/envs/pybind/lib/python3.9/site-packages/pybind11/include/pybind11/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(pybind11_DIR /Users/lg9783/miniconda3/envs/pybind/share/cmake/pybind11)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/sugarPos.csv
${PROJECT_SOURCE_DIR}/runClustering.py
${PROJECT_SOURCE_DIR}/test.csv
${PROJECT_SOURCE_DIR}/testFile.root
)

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

pybind11_add_module(clustering pyClustering.cc  ${sources} ${headers})
