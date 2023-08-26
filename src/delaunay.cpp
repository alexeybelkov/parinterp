//
// Created by lbelkov on 8/11/23.
//

#include "delaunay.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <parlay/sequence.h>
#include <Python.h>
#include "triangulation.h"

namespace py = pybind11;
using namespace delaunator;

Triangulation numpy_delaunay(const py::array_t<double, py::array::c_style | py::array::forcecast>& array) {
    uint32_t n = array.shape()[0];
    std::vector<double> P(2 * n);
    parlay::parallel_for(0, n, [&](uint32_t i) {
       P[2 * i] = array.at(i, 0);
       P[2 * i + 1] = array.at(i, 1);
    });
    Delaunator D(P);
    return {D};
}

#include "interpolator.h"

PYBIND11_MODULE(pydelaunay, m) {
    m.def("delaunay", &numpy_delaunay);
    py::class_<Triangulation>(m, "Triangulation", "Triangulation class")
            .def(py::init<>(), "Default constructor, does nothing")
            .def("find_triangle_naive", &Triangulation::find_triangle_naive,
                 "Return index of triangle, containing point p, else return -1", py::arg("x"), py::arg("y"))
            .def("find_triangle_bruteforce", &Triangulation::find_triangle_bruteforce,
                 "Return index of triangle, containing point p, else return -1", py::arg("x"), py::arg("y"))
            .def("jump_and_walk", &Triangulation::find_triangle_jump_and_walk,
                 "Using nearest vertex return index of triangle, containing point p, else return -1", py::arg("x"), py::arg("y"), py::arg("neighbor"))
            .def_readonly("vertices", &Triangulation::vertices, "numpy array of vertices")
            .def_readonly("triangles", &Triangulation::triangles, "numpy array of triangles");
    py::class_<BiLinearInterpolator>(m, "BiLinearInterpolator", "BiLinearInterpolator class")
            .def(py::init<const py::array_t<int32_t, py::array::c_style | py::array::forcecast>&,
                    const py::array_t<float, py::array::c_style | py::array::forcecast>&>())
            .def("__call__", &BiLinearInterpolator::call)
            .def_readonly("triangulation", &BiLinearInterpolator::triangulation)
//            .def("print_bar_coords", &BiLinearInterpolator::print_bar_coords)
            .def_readonly("values", &BiLinearInterpolator::values);
}

int main() {}