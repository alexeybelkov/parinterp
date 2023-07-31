#ifndef PARDELAUNAY_TRIANGULATION_H
#define PARDELAUNAY_TRIANGULATION_H

#endif //PARDELAUNAY_TRIANGULATION_H

class Triangulation {
private:
    pbbsbench::triangles<pbbsbench::pointT> triangulation;
public:
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> vertices;
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> triangles;
    Triangulation() = default;
    Triangulation(pbbsbench::triangles<pbbsbench::pointT>& triangulation);
    ~Triangulation();
    bool point_in_triangle(pbbsbench::pointT& p, pbbsbench::tri& triangle);
    bool point_in_triangle(int32_t& x, int32_t& y, pbbsbench::tri& triangle);
    int32_t find_triangle(pbbsbench::pointT& p);
    int32_t find_triangle(int32_t& x, int32_t& y);

};

Triangulation::Triangulation(pbbsbench::triangles<pbbsbench::pointT>& triangulation) {
    this -> triangulation = triangulation;
    uint32_t n = triangulation.P.size();
    uint32_t m = triangulation.T.size();
    std::vector<int64_t> points_shape = {long(triangulation.P.size()), 2};
    std::vector<int64_t> triangles_shape = {long(triangulation.T.size()), 3};
    std::vector<int64_t> points_strides = {sizeof(int32_t) * 2, sizeof(int32_t)};
    std::vector<int64_t> triangles_strides = {sizeof(int32_t) * 3, sizeof(int32_t)};
    std::vector<int32_t> flatten_points(2 * n);
    std::vector<int32_t> flatten_triangles(3 * m);
    parlay::parallel_for(0, n, [&](uint32_t i) {
        flatten_points[2 * i] = triangulation.P[i].x;
        flatten_points[2 * i + 1] = triangulation.P[i].y;
        });
    parlay::parallel_for(0, m, [&](uint32_t i) {
        for (uint32_t j = 0; j < 3; ++j)
            flatten_triangles[3 * i + j] = triangulation.T[i][j];
    });
    this->vertices = {std::move(points_shape), std::move(points_strides), flatten_points.data()};
    this->triangles = {std::move(triangles_shape), std::move(triangles_strides), flatten_triangles.data()};
}

Triangulation::~Triangulation() {
    this -> triangulation.P.clear(); // ???
}

bool Triangulation::point_in_triangle(pbbsbench::pointT& p, pbbsbench::tri& triangle) {
    pbbsbench::pointT A = triangulation.P[triangle[0]];
    pbbsbench::pointT B = triangulation.P[triangle[1]];
    pbbsbench::pointT C = triangulation.P[triangle[2]];
    auto s = (A.x - C.x) * (p.x - C.y) - (A.y - C.y) * (p.x - C.x);
    auto t = (B.x - A.x) * (p.y - A.y) - (B.y - A.y) * (p.x - A.x);

    if ((s < 0) != (t < 0) && s != 0 && t != 0) {
        return false;
    }

    auto d = (C.x - B.x) * (p.y - B.y) - (C.y - B.y) * (p.x - B.x);
    return d == 0 || (d < 0) == (s + t <= 0);
}

bool Triangulation::point_in_triangle(int32_t& x, int32_t& y, pbbsbench::tri& triangle) {
    pbbsbench::pointT A = triangulation.P[triangle[0]];
    pbbsbench::pointT B = triangulation.P[triangle[1]];
    pbbsbench::pointT C = triangulation.P[triangle[2]];
    auto s = (A.x - C.x) * (x - C.y) - (A.y - C.y) * (x - C.x);
    auto t = (B.x - A.x) * (y - A.y) - (B.y - A.y) * (x - A.x);

    if ((s < 0) != (t < 0) && s != 0 && t != 0) {
        return false;
    }

    auto d = (C.x - B.x) * (y - B.y) - (C.y - B.y) * (x - B.x);
    return d == 0 || (d < 0) == (s + t <= 0);
}

int32_t Triangulation::find_triangle(pbbsbench::pointT& p) {
    for (int32_t i = 0; i < triangulation.T.size(); ++i) {
        if (point_in_triangle(p, triangulation.T[i])) {
            return i;
        }
    }
    return -1;
}

int32_t Triangulation::find_triangle(int32_t& x, int32_t& y) {
    for (int32_t i = 0; i < triangulation.T.size(); ++i) {
        if (point_in_triangle(x, y, triangulation.T[i])) {
            return i;
        }
    }
    return -1;
}
