#ifndef PARDELAUNAY_TRIANGULATION_H
#define PARDELAUNAY_TRIANGULATION_H

#endif //PARDELAUNAY_TRIANGULATION_H

class Triangulation {
private:
    pbbsbench::triangles<pbbsbench::pointT> triangulation;
public:
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> vertices;
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> triangles;

    Triangulation(pbbsbench::triangles<pbbsbench::pointT> triangulation);
    uint32_t find_simplex(pbbsbench::pointT p);

};

Triangulation::Triangulation(pbbsbench::triangles<pbbsbench::pointT> triangulation) {
    this -> triangulation = triangulation;
    uint32_t n = triangulation.P.size();
    std::vector<long> points_shape = {long(triangulation.P.size()), 2};
    std::vector<long> triangles_shape = {long(triangulation.P.size()), 3};
    std::vector<long> points_strides = {sizeof(int32_t) * 2, sizeof(int32_t)};
    std::vector<long> triangles_strides = {sizeof(int32_t) * 3, sizeof(int32_t)};
    std::vector<int32_t> flatten_points(2 * n);
    std::vector<int32_t> flatten_triangles(3 * n);
    parlay::parallel_for(0, n, [&](uint32_t i) {
        flatten_points[2 * i] = int32_t(triangulation.P[i].x);
        flatten_points[2 * i + 1] = int32_t(triangulation.P[i].y);
        for (uint32_t j = 0; j < 3; ++j)
            flatten_triangles[3 * i + j] = triangulation.T[i][j];
        });
    this->vertices = {std::move(points_shape), std::move(points_strides), flatten_points.data()};
    this->triangles = {std::move(triangles_shape), std::move(triangles_strides), flatten_triangles.data()};
}

uint32_t Triangulation::find_simplex(pbbsbench::pointT p) {
    parlay::parallel_for(0, this->triangles.size(), [&] (uint32_t i) {
        this -> triangulation.T[i];
    });
}