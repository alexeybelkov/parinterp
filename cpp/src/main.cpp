#pragma GCC target("avx2")

#include <iostream>
#include <stdexcept>
#include "triangulator.h"
#include "interpolator.h"
#include <Python.h>

PYBIND11_MODULE(parinterp, m) {
    py::class_<Interpolator>(m, "Linear2DInterpolatorCpp", "Interpolator class")
        .def(py::init<const Triangulator::pyarr_size_t&, int>())
        .def("__call__", &Interpolator::operator())
        .def("naive_points_location", &Interpolator::naive_points_location);
}

int main() {
}
