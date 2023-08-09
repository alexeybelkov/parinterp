// The inteface for delaunay triangulation
#include "geometry.h"
#include "parlay/primitives.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include "triangulation.h"

class Triangulation;

namespace py = pybind11;

namespace pbbsbench {

  using coord = double;
  using pointT = point2d<coord>;

  triangles<pointT> delaunay(const parlay::sequence<pointT>& P);
  Triangulation numpy_delaunay(const py::array_t<double, py::array::c_style | py::array::forcecast>& array);
}