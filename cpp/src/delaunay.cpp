#include <iostream>
#include <stdexcept>
#include "triangulation.h"
#include "interpolator.h"
#include <Python.h>

// PYBIND11_MODULE(pydelaunay, m) {
//     py::class_<Interpolator>(Linear2DInterpolator)
// }


// #include "interpolator.h"

// PYBIND11_MODULE(pydelaunay, m) {
//     m.def("delaunay", &numpy_delaunay);
//     py::class_<Triangulation>(m, "Triangulation", "Triangulation class")
//             .def(py::init<>(), "Default constructor, does nothing")
//             .def("find_triangle_naive", &Triangulation::find_triangle_naive,
//                  "Return index of triangle, containing point p, else return -1", py::arg("x"), py::arg("y"))
//             .def("find_triangle_bruteforce", &Triangulation::find_triangle_bruteforce,
//                  "Return index of triangle, containing point p, else return -1", py::arg("x"), py::arg("y"))
//             .def("jump_and_walk", &Triangulation::find_triangle_jump_and_walk,
//                  "Using nearest vertex return index of triangle, containing point p, else return -1", py::arg("x"), py::arg("y"), py::arg("neighbor"))
//             .def_readonly("vertices", &Triangulation::vertices, "numpy array of vertices")
//             .def_readonly("triangles", &Triangulation::triangles, "numpy array of triangles");
//     py::class_<BiLinearInterpolator>(m, "BiLinearInterpolator", "BiLinearInterpolator class")
//             .def(py::init<const py::array_t<int32_t, py::array::c_style | py::array::forcecast>&,
//                     const py::array_t<float, py::array::c_style | py::array::forcecast>&>())
//             .def("__call__", &BiLinearInterpolator::call)
//             .def_readonly("triangulation", &BiLinearInterpolator::triangulation)
//             .def_readonly("values", &BiLinearInterpolator::values);
// }

py::array_t<size_t, py::array::c_style | py::array::forcecast> make_pyarray(size_t size) {
    // py::scoped_interpreter guard{};
    std::vector<size_t> v(size);
    // #pragma parallel for
    for (size_t i = 0; i < size; ++i) {
        v[i] = static_cast<int>(i * i);
    }
    std::vector<int64_t> shape = {int64_t(size / 2), int64_t(2)};
    std::vector<int64_t> strides = {int64_t(sizeof(size_t)), int64_t(sizeof(int))};
    return {std::move(shape), std::move(strides), v.data()};
}

int main() {
    // omp_set_num_threads(4);
    py::scoped_interpreter guard{};
    auto arr = make_pyarray(100);
    // delete[] &guard;
    Triangulator tri(arr, 1);
    std::cout << tri.triangles.size() << std::endl;
    std::cout << "TRIANGULATED\n";

    Interpolator inter(arr, 1);
    std::vector<double> vec_values(arr.shape()[0]);
    for (size_t i = 0; i < vec_values.size(); ++i) {
        vec_values[i] = static_cast<double>(i);
    }
    // Interpolator::pyarr_double values{{vec_values.size()}, {sizeof(double)}, vec_values.data()};
    // inter(values, arr, arr);
}