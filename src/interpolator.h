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
    float bilinear_barycentric_interpolation(int32_t& x, int32_t& y, int32_t& t);
    py::array_t<float, py::array::c_style | py::array::forcecast> call(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points);

};

BiLinearInterpolator::BiLinearInterpolator(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points,
                                           const py::array_t<float, py::array::c_style | py::array::forcecast>& values, float fill_value) {
    this -> values = values;
    this -> triangulation = pbbsbench::numpy_delaunay(points);
    this -> fill_value = fill_value;
}

float BiLinearInterpolator::bilinear_barycentric_interpolation(int32_t &x, int32_t &y, int32_t& t) {

    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Vertex_approach

    int32_t i = this -> triangulation.triangles.at(t, 0);
    int32_t j = this -> triangulation.triangles.at(t, 1);
    int32_t k = this -> triangulation.triangles.at(t, 2);

    int32_t x1 = triangulation.vertices.at(i, 0);
    int32_t y1 = triangulation.vertices.at(i, 1);
    int32_t x2 = triangulation.vertices.at(j, 0);
    int32_t y2 = triangulation.vertices.at(j, 1);
    int32_t x3 = triangulation.vertices.at(k, 0);
    int32_t y3 = triangulation.vertices.at(k, 1);

    int32_t dx32 = x3 - x2;
    int32_t dx13 = x1 - x3;
    int32_t dx21 = x2 - x1;

    int32_t dy23 = y2 - y3;
    int32_t dy31 = y3 - y1;
    int32_t dy12 = y1 - y2;

    int32_t dxy23 = x2*y3 - x3*y2;
    int32_t dxy31 = x3*y1 - x1*y3;
    int32_t dxy12 = x1*y2 - x2*y1;

    float lambda_1 = float(dxy23 + x * dy23 + y * dx32);
    float lambda_2 = float(dxy31 + x * dy31 + y * dx13);
    float lambda_3 = float(dxy12 + x * dy12 + y * dx21);

    float f1 = this -> values.at(i);
    float f2 = this -> values.at(j);
    float f3 = this -> values.at(k);

    return (lambda_1 * f1 + lambda_2 * f2 + lambda_3 * f3) / float(dxy12 + dxy31 + dxy23);
}

py::array_t<float, py::array::c_style | py::array::forcecast> BiLinearInterpolator::call(const py::array_t<int32_t, py::array::c_style | py::array::forcecast>& points) {
    uint32_t n = points.shape()[0];
    std::vector<float> interpolated(n, this -> fill_value);

    parlay::parallel_for(0, n, [&](uint32_t i) {
        int32_t x = points.at(i, 0);
        int32_t y = points.at(i, 1);
        int32_t t = this -> triangulation.find_triangle(x, y);
        if (t != -1)
            interpolated[i] = bilinear_barycentric_interpolation(x, y, t);
    });

    std::vector<int64_t> shape = {n, 2};
    std::vector<int64_t> strides = {sizeof(float), sizeof(float)};

    return {std::move(shape), std::move(strides), interpolated.data()};
}
