cmake_minimum_required(VERSION 3.23)
project(pardelaunay)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_BUILD_TYPE Debug)

# find_package(OpenMP)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
# endif()

# include_directories(/homes/lbelkov/packages/cmake/parlaylib/include)
# include_directories(/include/)
include_directories(/homes/lbelkov/miniconda3/envs/py3104/include/python3.10)
include_directories(/homes/lbelkov/miniconda3/envs/py3104/include/)
include_directories(/homes/lbelkov/miniconda3/envs/py3104/lib/python3.10/site-packages/pybind11/include)
set(PYTHON_LIBRARY /homes/lbelkov/miniconda3/envs/py3104/lib/libpython3.10.so)

set(pybind11_DIR /homes/lbelkov/miniconda3/envs/py3104/lib/python3.10/site-packages/pybind11/share/cmake/pybind11)
# set(PARLAY_NUM_THREADS 16)


find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(pydelaunay src/delaunay.cpp)
add_executable(pardelaunay src/delaunay.cpp)

target_link_libraries(pardelaunay ${PYTHON_LIBRARIES})
target_link_libraries(pardelaunay ${pybind11_LIBRARIES})
