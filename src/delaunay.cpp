//
// Created by lbelkov on 8/5/23.
//

#include "Python.h"
#include "delaunay.h"
#include "triangulation.h"
#include "interpolator.h"
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;


py::array_t<int32_t, py::array::c_style | py::array::forcecast> numpy_delaunay(const py::array_t<double, py::array::c_style | py::array::forcecast>& array) {
    uint32_t n = array.shape()[0];
    auto P = parlay::tabulate(n, [&](int32_t i) -> point {
        return {i, array.at(i, 0), array.at(i, 1)};
    });
    Delaunay DT(P);
    auto triangles = DT.mesh.keys();
    parlay::sequence<int32_t> flatten_triangles(3 * triangles.size());
    parlay::parallel_for(0, n, [&] (uint32_t i) {
       flatten_triangles[3 * i] = triangles[i][0];
       flatten_triangles[3 * i + 1] = triangles[i][1];
       flatten_triangles[3 * i + 2] = triangles[i][2];
    });
    std::vector<uint32_t> strides = {sizeof(int32_t) * 3, sizeof(int32_t)};
    std::vector<uint32_t> shape = {uint32_t(triangles.size()), 3};
    return {std::move(shape), std::move(strides), flatten_triangles.data()};
}


PYBIND11_MODULE(pydelaunay, m) {
    m.def("delaunay", &numpy_delaunay);
}

int main() {}
