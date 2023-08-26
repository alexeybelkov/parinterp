//
// Created by lbelkov on 7/31/23.
//

#ifndef PARDELAUNAY_INTERPOLATOR_H
#define PARDELAUNAY_INTERPOLATOR_H

#endif //PARDELAUNAY_INTERPOLATOR_H

class BiLinearInterpolator {

public:
    Triangulation triangulation;
    py::array_t<float, py::array::c_style | py::array::forcecast> values;
//    py::array_t<float, py::array::c_style | py::array::forcecast> barycenters;
    BiLinearInterpolator() = default;
    BiLinearInterpolator(const py::array_t<double, py::array::c_style | py::array::forcecast>& points,
                         const py::array_t<float, py::array::c_style | py::array::forcecast>& values);
    double bilinear_barycentric_interpolation(const int64_t t, const parlay::sequence<double>& bar_coords);
    py::array_t<float, py::array::c_style | py::array::forcecast> call(const py::array_t<double, py::array::c_style | py::array::forcecast>& points,
                                                                       const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& neighbors,
                                                                       float fill_value = 0.0);
//    void print_bar_coords(double x, double y, uint32_t t) {
//        auto bar_coords = this -> triangulation.barycentric_coordinates(x, y, t);
//        for (auto& b: bar_coords) {
//            auto r1 = b.first;
//            auto i = b.second;
//            auto f1 = i < triangulation.triangulation.numPoints() - triangulation.boundary_size ? values.at(i) : 0.0;
//            std::cout << r1 << ' ' << i << ' ' << f1 << '\n';
//        }
//    }
};

BiLinearInterpolator::BiLinearInterpolator(const py::array_t<double, py::array::c_style | py::array::forcecast>& points,
                                           const py::array_t<float, py::array::c_style | py::array::forcecast>& values) {
    this -> values = values;
    this -> triangulation = numpy_delaunay(points);
}

double BiLinearInterpolator::bilinear_barycentric_interpolation(const int64_t t, const parlay::sequence<double>& bar_coords) {
    uint32_t i = triangulation.triangles.at(t, 0);
    uint32_t j = triangulation.triangles.at(t, 1);
    uint32_t k = triangulation.triangles.at(t, 2);
    double f1 = values.at(i);
    double f2 = values.at(j);
    double f3 = values.at(k);
    auto interpolated = (bar_coords[0] * f1 + bar_coords[1] * f2 + bar_coords[2] * f3);
    return interpolated;
}

py::array_t<float, py::array::c_style | py::array::forcecast> BiLinearInterpolator::call(const py::array_t<double, py::array::c_style | py::array::forcecast>& points,
                                                                                         const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& neighbors,
                                                                                         float fill_value) {
    uint32_t n = points.shape()[0];
    parlay::sequence<float> interpolated(n, fill_value);

    parlay::parallel_for(0, n, [&](uint32_t i) {
        double x = points.at(i, 0);
        double y = points.at(i, 1);

//        int32_t t = triangulation.find_triangle_bruteforce(x, y);
//        if (t != -1) {
//            auto bar_coords = triangulation.barycentric_coordinates(x, y, uint32_t(t));
//            interpolated[i] = bilinear_barycentric_interpolation(t, bar_coords);
//        }

       auto neighbor = neighbors.at(i);
       auto checked = triangulation.jump_and_walk(x, y, neighbor);
       if (checked.first != -1) {
           interpolated[i] = bilinear_barycentric_interpolation(checked.first, checked.second);
       }
    });

    std::vector<int64_t> shape = {n};
    std::vector<int64_t> strides = {sizeof(float)};

    return {std::move(shape), std::move(strides), interpolated.data()};
}
