//
// Created by lbelkov on 8/11/23.
//

#include "delaunay.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <parlay/sequence.h>
#include <Python.h>

namespace py = pybind11;

class Triangulation {
public:
    py::array_t<double, py::array::c_style | py::array::forcecast> vertices;
    py::array_t<size_t, py::array::c_style | py::array::forcecast> triangles;
    Triangulation() = default;
    Triangulation(delaunator::Delaunator &D) {
        uint32_t n = D.coords.size();
        uint32_t m = D.triangles.size();
        std::vector<int64_t> points_shape = {n / 2, 2};
        std::vector<int64_t> triangles_shape = {m / 3, 3};
        std::vector<int64_t> points_strides = {sizeof(double) * 2, sizeof(double)};
        std::vector<int64_t> triangles_strides = {sizeof(size_t) * 3, sizeof(size_t)};

        this->vertices = {std::move(points_shape), std::move(points_strides), D.coords.data()};
        this->triangles = {std::move(triangles_shape), std::move(triangles_strides), D.triangles.data()};
    }
};

Triangulation numpy_delaunay(const py::array_t<double, py::array::c_style | py::array::forcecast>& array) {
    uint32_t n = array.shape()[0];
    std::vector<double> P(2 * n);
    parlay::parallel_for(0, n, [&](uint32_t i) {
       P[2 * i] = array.at(i, 0);
       P[2 * i + 1] = array.at(i, 1);
    });
    delaunator::Delaunator D(P);
    return {D};
}


PYBIND11_MODULE(pydelaunay, m) {
    m.def("delaunay", &numpy_delaunay);
    py::class_<Triangulation>(m, "Triangulation", "Triangulation class")
            .def(py::init<>(), "Default constructor, does nothing")
            .def_readonly("vertices", &Triangulation::vertices, "numpy array of vertices")
            .def_readonly("triangles", &Triangulation::triangles, "numpy array of triangles");
}

int main() {}