#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <omp.h>

int main() {
    size_t NUM_THREADS;
    std::cin >> NUM_THREADS;
    std::cout << omp_get_num_threads() << std::endl;
    omp_set_num_threads(NUM_THREADS);
    std::cout << omp_get_num_threads() << std::endl;
    #pragma omp parallel
    {
        std::cout << "Hello\n";
    }
    return 0;
}