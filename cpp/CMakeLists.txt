cmake_minimum_required(VERSION 3.23)
project(parinterp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fopenmp")
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)
set(CMAKE_BUILD_TYPE Release)

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})

# find_package(OpenMP)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
# endif()

# include_directories(/homes/lbelkov/packages/cmake/parlaylib/include)
# include_directories(/include/)
# include_directories(/homes/lbelkov/miniconda3/envs/py3104/include/python3.10)
# include_directories(/homes/lbelkov/miniconda3/envs/py3104/include/)
# include_directories(/homes/lbelkov/miniconda3/envs/py3104/lib/python3.10/site-packages/pybind11/include)
# set(PYTHON_LIBRARY /homes/lbelkov/miniconda3/envs/py3104/lib/libpython3.10.so)

# set(pybind11_DIR /homes/lbelkov/miniconda3/envs/py3104/lib/python3.10/site-packages/pybind11/share/cmake/pybind11)
# set(PARLAY_NUM_THREADS 16)

# include_directories(${PYTHON_LIBRARIE})

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(parinterp src/main.cpp)
# target_link_libraries(parinterp PUBLIC)
# add_executable(main src/main.cpp)
# target_link_libraries(main ${PYTHON_LIBRARIES})
# target_link_libraries(main ${pybind11_LIBRARIES})

# install(TARGETS parinterp COMPONENT python LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}")
