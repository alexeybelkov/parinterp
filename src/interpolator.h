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
    float fill_value = 0.0;
    BiLinearInterpolator() = default;
    BiLinearInterpolator(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points,
                         const py::array_t<float, py::array::c_style | py::array::forcecast>& values, float fill_value = 0.0);
    float interpolate(int32_t& x, int32_t& y, int32_t& t);
    py::array_t<float, py::array::c_style | py::array::forcecast> call(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points);

};

BiLinearInterpolator::BiLinearInterpolator(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points,
                                           const py::array_t<float, py::array::c_style | py::array::forcecast>& values, float fill_value) {
    this -> values = values;
    this -> triangulation = pbbsbench::numpy_delaunay(points);
    this -> fill_value = fill_value;
}

float BiLinearInterpolator::interpolate(int32_t &x, int32_t &y, int32_t& t) {
    int32_t i = this -> triangulation.triangles.at(t, 0);
    int32_t j = this -> triangulation.triangles.at(t, 1);
    int32_t k = this -> triangulation.triangles.at(t, 2);
    return this -> fill_value;
}

py::array_t<float, py::array::c_style | py::array::forcecast> BiLinearInterpolator::call(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points) {
    uint32_t n = points.shape()[0];
    std::vector<float> interpolated(n, this -> fill_value);

    parlay::parallel_for(0, n, [&](uint32_t i) {
        int32_t x = points.at(i, 0);
        int32_t y = points.at(i, 1);
        int32_t t = this -> triangulation.find_triangle(x, y);
        if (t != -1)
            interpolated[i] = interpolate(x, y, t);
    });

    std::vector<int64_t> shape = {n, 2};
    std::vector<int64_t> strides = {sizeof(int32_t), sizeof(int32_t)};

    return {std::move(shape), std::move(strides), interpolated.data()};
}
