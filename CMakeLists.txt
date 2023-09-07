cmake_minimum_required(VERSION 3.25)
project(pardelaunay)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(/homes/lbelkov/packages/cmake/parlaylib/include)
#include_directories(/include/)
include_directories(/homes/lbelkov/miniconda3/envs/py3104/include/python3.10)
include_directories(/homes/lbelkov/miniconda3/envs/py3104/include/)
include_directories(/homes/lbelkov/miniconda3/envs/py3104/lib/python3.10/site-packages/pybind11/include)
set(PYTHON_LIBRARY /homes/lbelkov/miniconda3/envs/py3104/lib/libpython3.10.so)

set(pybind11_DIR /homes/lbelkov/miniconda3/envs/py3104/lib/python3.10/site-packages/pybind11/share/cmake/pybind11)
set(PARLAY_NUM_THREADS 16)
###


find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(pydelaunay src/delaunay.cpp)
add_executable(pardelaunay src/delaunay.cpp)
#        src/interpolator.h)
#        pargeo/as_numpy.h)

target_link_libraries(pardelaunay ${PYTHON_LIBRARIES})
