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
    py::array_t<float, py::array::c_style | py::array::forcecast> barycenters;
    BiLinearInterpolator() = default;
    BiLinearInterpolator(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points,
                         const py::array_t<float, py::array::c_style | py::array::forcecast>& values);
    double bilinear_barycentric_interpolation(int32_t& x, int32_t& y, int32_t& t);
    py::array_t<float, py::array::c_style | py::array::forcecast> call(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points,
                                                                       const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& neighbors,
                                                                       float fill_value = 0.0);
};

BiLinearInterpolator::BiLinearInterpolator(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points,
                                           const py::array_t<float, py::array::c_style | py::array::forcecast>& values) {
    this -> values = values;
    this -> triangulation = pbbsbench::numpy_delaunay(points);

    uint32_t m = this -> triangulation.triangulation.T.size();

    parlay::sequence<float> barycenters_(2 * m);

    parlay::parallel_for(0, m, [&](uint32_t t) {
        int32_t i = this -> triangulation.triangulation.T[t][0];
        int32_t j = this -> triangulation.triangulation.T[t][1];
        int32_t k = this -> triangulation.triangulation.T[t][2];

        int32_t x1 = triangulation.triangulation.P[i].x;
        int32_t y1 = triangulation.triangulation.P[i].y;
        int32_t x2 = triangulation.triangulation.P[j].x;
        int32_t y2 = triangulation.triangulation.P[j].y;
        int32_t x3 = triangulation.triangulation.P[k].x;
        int32_t y3 = triangulation.triangulation.P[k].y;

        float x = float(x1 + x2 + x3) / 3.0;
        float y = float(y1 + y2 + y3) / 3.0;
        barycenters_[2 * t] = x;
        barycenters_[2 * t + 1] = y;
    });

    std::vector<int64_t> shape = {m, 2};
    std::vector<int64_t> strides = {sizeof(float) * 2, sizeof(float)};
    this -> barycenters = {std::move(shape), std::move(strides), barycenters_.data()};

}

double BiLinearInterpolator::bilinear_barycentric_interpolation(int32_t &x, int32_t &y, int32_t& t) {

    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Vertex_approach

    int32_t i = triangulation.triangulation.T[t][0];
    int32_t j = triangulation.triangulation.T[t][1];
    int32_t k = triangulation.triangulation.T[t][2];

    int32_t x1 = triangulation.triangulation.P[i].x;
    int32_t y1 = triangulation.triangulation.P[i].y;
    int32_t x2 = triangulation.triangulation.P[j].x;
    int32_t y2 = triangulation.triangulation.P[j].y;
    int32_t x3 = triangulation.triangulation.P[k].x;
    int32_t y3 = triangulation.triangulation.P[k].y;

    int32_t dx32 = x3 - x2;
    int32_t dx13 = x1 - x3;
    int32_t dx21 = x2 - x1;

    int32_t dy23 = y2 - y3;
    int32_t dy31 = y3 - y1;
    int32_t dy12 = y1 - y2;

    int32_t dxy23 = x2*y3 - x3*y2;
    int32_t dxy31 = x3*y1 - x1*y3;
    int32_t dxy12 = x1*y2 - x2*y1;

    auto lambda_1 = double(dxy23 + x * dy23 + y * dx32);
    auto lambda_2 = double(dxy31 + x * dy31 + y * dx13);
    auto lambda_3 = double(dxy12 + x * dy12 + y * dx21);

    double f1 = this -> values.at(i);
    double f2 = this -> values.at(j);
    double f3 = this -> values.at(k);

    auto interpolated = (lambda_1 * f1 + lambda_2 * f2 + lambda_3 * f3) / double(dxy12 + dxy31 + dxy23);
    return interpolated;
}

py::array_t<float, py::array::c_style | py::array::forcecast> BiLinearInterpolator::call(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points,
                                                                                         const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& neighbors,
                                                                                         float fill_value) {
    uint32_t n = points.shape()[0];
    parlay::sequence<float> interpolated(n, fill_value);

    parlay::parallel_for(0, n, [&](uint32_t i) {
        int32_t x = points.at(i, 0);
        int32_t y = points.at(i, 1);
        int32_t neighbor = neighbors.at(i);
        int32_t t = this -> triangulation.check_neighborhood(x, y, neighbor);
        if (t != -1)
            interpolated[i] = bilinear_barycentric_interpolation(x, y, t);
    });

    std::vector<int64_t> shape = {n};
    std::vector<int64_t> strides = {sizeof(float)};

    return {std::move(shape), std::move(strides), interpolated.data()};
}
