// The inteface for delaunay triangulation
#include "geometry.h"
#include "parlay/primitives.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace pbbsbench {

  using coord = double;
  using pointT = point2d<coord>;

  triangles<pointT> delaunay(parlay::sequence<pointT>& P);
  py::array_t<int32_t, py::array::c_style | py::array::forcecast> numpy_delaunay(py::array_t<int32_t, py::array::c_style | py::array::forcecast>& array);
}