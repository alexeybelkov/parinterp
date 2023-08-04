#ifndef PARDELAUNAY_TRIANGULATION_H
#define PARDELAUNAY_TRIANGULATION_H

#endif //PARDELAUNAY_TRIANGULATION_H

uint32_t bijection(uint32_t& a, uint32_t& b) {
    return (a + b) * (a + b + 1) / 2 + b; // Bijection: NxN -> N
}

class Triangulation {
public:
//    parlay::sequence<std::pair<int32_t, int32_t>> edges;
    pbbsbench::triangles<pbbsbench::pointT> triangulation;
    parlay::sequence<parlay::sequence<uint32_t>> adj_p2t;
    std::unordered_map<uint32_t, parlay::sequence<uint32_t>> adj_e2t;
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> vertices;
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> triangles;
    Triangulation() = default;
    Triangulation(pbbsbench::triangles<pbbsbench::pointT>& triangulation);
    ~Triangulation();
    bool point_in_triangle(parlay::sequence<pair<double, uint32_t>>& bar_coords);
    bool point_in_triangle(int32_t x, int32_t y, pbbsbench::tri& triangle);
    int32_t check_adj_triangles(int32_t x, int32_t y, uint32_t t);
    int32_t find_triangle(pbbsbench::pointT& p);
    int32_t find_triangle(int32_t x, int32_t y);
    int32_t check_neighborhood(int32_t x, int32_t y, int32_t pi);
    parlay::sequence<std::pair<double, uint32_t>> barycentric_coordinates(int32_t x, int32_t y, uint32_t t);
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

    this -> adj_p2t = parlay::sequence<parlay::sequence<uint32_t>>(n);

    for (uint32_t i = 0; i < m; ++i) {
        for (int32_t j = 0; j < 3; ++j) {
            uint32_t a = triangulation.T[i][j];
            uint32_t b = triangulation.T[i][(j + 1) % 3];

            adj_p2t[a].push_back(i);

            if (a > b)
                std::swap(a, b);
            uint32_t e = bijection(a, b);
            if (adj_e2t.find(e) == adj_e2t.end()) {
                parlay::sequence<uint32_t> s(1, i);
                adj_e2t.insert({e, s});
            }
            else {
                adj_e2t[e].push_back(i);
            }
        }
    }


    this->vertices = {std::move(points_shape), std::move(points_strides), flatten_points.data()};
    this->triangles = {std::move(triangles_shape), std::move(triangles_strides), flatten_triangles.data()};
}

Triangulation::~Triangulation() {
    this -> triangulation.P.clear(); // ???
    this -> adj_p2t.clear();
    this -> adj_e2t.clear();
}

bool Triangulation::point_in_triangle(parlay::sequence<pair<double, uint32_t>>& bar_coords) {
    return bar_coords[0].first > 0 and bar_coords[1].first > 0 and bar_coords[2].first > 0;
}

bool Triangulation::point_in_triangle(int32_t x, int32_t y, pbbsbench::tri& triangle) {
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

//int32_t Triangulation::find_triangle(pbbsbench::pointT& p) {
//    for (int32_t i = 0; i < triangulation.T.size(); ++i) {
//        if (point_in_triangle(p, triangulation.T[i])) {
//            return i;
//        }
//    }
//    return -1;
//}

int32_t Triangulation::find_triangle(int32_t x, int32_t y) {
    for (int32_t i = 0; i < triangulation.T.size(); ++i) {
        if (point_in_triangle(x, y, triangulation.T[i])) {
            return i;
        }
    }
    return -1;
}

int32_t Triangulation::check_adj_triangles(int32_t x, int32_t y, uint32_t t) {
    for (uint32_t i = 0; i < 3; ++i) {
        uint32_t a = triangulation.T[t][i];
        uint32_t b = triangulation.T[t][(i + 1) % 3];
        if (a > b)
            std::swap(a, b);

        uint32_t e = bijection(a, b);
        for (auto& tt : adj_e2t[e]) {
            if (tt != t) {
                auto bar_coords = barycentric_coordinates(x, y, tt);
                if (point_in_triangle(bar_coords))
                    return tt;
            }
        }
    }
    return -1;
}

int32_t Triangulation::check_neighborhood(int32_t x, int32_t y, int32_t pi) {
    for (auto& t : this -> adj_p2t[pi]) {   // A bit of overhead :retard:
        auto j = check_adj_triangles(x, y, t);
        if (j != -1)
            return j;
    }
    return -1;
}

parlay::sequence<std::pair<double, uint32_t>> Triangulation::barycentric_coordinates(int32_t x, int32_t y, uint32_t t) {
    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Vertex_approach

    uint32_t i = triangulation.T[t][0];
    uint32_t j = triangulation.T[t][1];
    uint32_t k = triangulation.T[t][2];

    int32_t x1 = triangulation.P[i].x;
    int32_t y1 = triangulation.P[i].y;
    int32_t x2 = triangulation.P[j].x;
    int32_t y2 = triangulation.P[j].y;
    int32_t x3 = triangulation.P[k].x;
    int32_t y3 = triangulation.P[k].y;

    int32_t dx32 = x3 - x2;
    int32_t dx13 = x1 - x3;
    int32_t dx21 = x2 - x1;

    int32_t dy23 = y2 - y3;
    int32_t dy31 = y3 - y1;
    int32_t dy12 = y1 - y2;

    int32_t dxy23 = x2*y3 - x3*y2;
    int32_t dxy31 = x3*y1 - x1*y3;
    int32_t dxy12 = x1*y2 - x2*y1;

    double one_over_2area = 1.0 / double(dxy12 + dxy31 + dxy23);

    auto lambda_1 = double(dxy23 + x * dy23 + y * dx32) * one_over_2area;
    auto lambda_2 = double(dxy31 + x * dy31 + y * dx13) * one_over_2area;
    auto lambda_3 = double(dxy12 + x * dy12 + y * dx21) * one_over_2area;
    return {{lambda_1, i}, {lambda_2, j}, {lambda_3, k}};
}
