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
    double bilinear_barycentric_interpolation(int32_t t, double fill_value, parlay::sequence<std::pair<double, uint32_t>>& bar_coords);
    py::array_t<float, py::array::c_style | py::array::forcecast> call(const py::array_t<double, py::array::c_style | py::array::forcecast>& points,
                                                                       const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& neighbors,
                                                                       float fill_value = 0.0);
    void print_bar_coords(double x, double y, uint32_t t) {
        auto bar_coords = this -> triangulation.barycentric_coordinates(x, y, t);
        for (auto& b: bar_coords) {
            auto r1 = b.first;
            auto i = b.second;
            auto f1 = i < triangulation.triangulation.numPoints() - triangulation.boundary_size ? values.at(i) : 0.0;
            std::cout << r1 << ' ' << i << ' ' << f1 << '\n';
        }
    }
};

BiLinearInterpolator::BiLinearInterpolator(const py::array_t<double, py::array::c_style | py::array::forcecast>& points,
                                           const py::array_t<float, py::array::c_style | py::array::forcecast>& values) {
    this -> values = values;
    this -> triangulation = pbbsbench::numpy_delaunay(points);

    uint32_t m = this -> triangulation.triangulation.T.size();

//    parlay::sequence<float> barycenters_(2 * m);
//
//    parlay::parallel_for(0, m, [&](uint32_t t) {
//        int32_t i = this -> triangulation.triangulation.T[t][0];
//        int32_t j = this -> triangulation.triangulation.T[t][1];
//        int32_t k = this -> triangulation.triangulation.T[t][2];
//
//        int32_t x1 = triangulation.triangulation.P[i].x;
//        int32_t y1 = triangulation.triangulation.P[i].y;
//        int32_t x2 = triangulation.triangulation.P[j].x;
//        int32_t y2 = triangulation.triangulation.P[j].y;
//        int32_t x3 = triangulation.triangulation.P[k].x;
//        int32_t y3 = triangulation.triangulation.P[k].y;
//
//        float x = float(x1 + x2 + x3) / 3.0;
//        float y = float(y1 + y2 + y3) / 3.0;
//        barycenters_[2 * t] = x;
//        barycenters_[2 * t + 1] = y;
//    });

//    std::vector<int64_t> shape = {m, 2};
//    std::vector<int64_t> strides = {sizeof(float) * 2, sizeof(float)};
//    this -> barycenters = {std::move(shape), std::move(strides), barycenters_.data()};

}

double BiLinearInterpolator::bilinear_barycentric_interpolation(int32_t t, double fill_value, parlay::sequence<std::pair<double, uint32_t>>& bar_coords) {
    uint32_t i = bar_coords[0].second;
    uint32_t j = bar_coords[0].second;
    uint32_t k = bar_coords[0].second;
    uint32_t n = this -> triangulation.triangulation.numPoints() - this -> triangulation.boundary_size;
    double f1 = i < n ? this -> values.at(i) : fill_value;
    double f2 = j < n ? this -> values.at(j) : fill_value;
    double f3 = k < n ? this -> values.at(k) : fill_value;
    auto interpolated = (bar_coords[0].first * f1 + bar_coords[1].first * f2 + bar_coords[2].first * f3);
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
        int32_t neighbor = neighbors.at(i);
        int32_t t = this -> triangulation.check_neighborhood(x, y, neighbor);
        if (t != -1) {
            auto bar_coords = triangulation.barycentric_coordinates(x, y, t);
            interpolated[i] = bilinear_barycentric_interpolation(t, fill_value, bar_coords);
        }
    });

    std::vector<int64_t> shape = {n};
    std::vector<int64_t> strides = {sizeof(float)};

    return {std::move(shape), std::move(strides), interpolated.data()};
}
