#ifndef PARDELAUNAY_TRIANGULATION_H
#define PARDELAUNAY_TRIANGULATION_H

#endif //PARDELAUNAY_TRIANGULATION_H

#include <parlay/primitives.h>

// Szudzik’s Pairing Function https://en.wikipedia.org/wiki/Pairing_function#Other_pairing_functions
uint64_t elegant_pair(const int32_t& a, const int32_t& b) {
    auto aa = uint64_t(a);
    auto bb = uint64_t(b);
    return bb * bb + aa;
}

uint64_t elegant_pair(const uint32_t& a, const uint32_t& b) {
    auto aa = uint64_t(a);
    auto bb = uint64_t(b);

    return bb * bb + aa;
//    return (aa << 15) + bb;
//    if (aa > bb)
//        return bb*bb + aa;
//    else
//        return aa * aa + aa + bb;
}

uint64_t elegant_pair(const std::pair<int32_t, int32_t>& e) {
    uint64_t aa = e.first;
    uint64_t bb = e.second;
    return bb * bb + aa;
//    return (aa << 15) + bb;
//    auto p = (a + b);
//    return (p * (p + 1) + 2 * a) / 2;
//    if (aa > bb)
//        return bb*bb + aa;
//    else
//        return aa * aa + aa + bb;
}

int32_t sign(const int32_t& x) {
    return x >= 0 ? x ? 1 : 0 : -1;
}


bool edges_eq(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    return (a == c) and (b == d);
}

// https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line_segment

bool segm_intersection(double Ax, double Ay, double Bx, double By,
                       double Cx, double Cy, double Dx, double Dy) {

    auto AxCy = Ax * Cy;
    auto AxDy = Ax * Dy;
    auto CxAy = Cx * Ay;
    auto DxAy = Dx * Ay;
    auto CxBy = Cx * By;
    auto BxCy = Bx * Cy;

//    auto unnorm_t = AxCy - AxDy + Cx * Dy - CxAy + DxAy - Dx * Cy;
//    auto unnorm_u = -Ax * By - CxAy + CxBy + Bx * Ay + AxCy - BxCy;
//    auto denum = AxCy - AxDy - BxCy + Bx * Dy - CxAy + DxAy - Dx * By + CxBy;

    auto unnorm_t = (Ax - Cx) * (Cy - Dy) - (Ay - Cy) * (Cx - Dx);
    auto unnorm_u = (Ax - Cx) * (Ay - By) - (Ay - Cy) * (Ax - Bx);
    auto denum = (Ax - Bx) * (Cy - Dy) - (Ay - By) * (Cx - Dx);

    if (denum == 0.0)
        return false;

    bool is_t = unnorm_t >= 0 and unnorm_t <= denum;
    bool is_u = unnorm_u >= 0 and unnorm_u <= denum;
    return is_t and is_u;
}

namespace py = pybind11;
using namespace delaunator;

class Triangulation {
public:
//    Delaunator* delauany;
    uint32_t numVertices;
    uint32_t numTriangles;
    parlay::sequence<parlay::sequence<uint32_t>> adj_p2t;
    std::unordered_map<uint64_t, parlay::sequence<uint32_t>> adj_e2t;
//    std::vector<std::vector<std::vector<uint32_t>>> adj_e2t;
    py::array_t<double, py::array::c_style | py::array::forcecast> vertices;
    py::array_t<size_t, py::array::c_style | py::array::forcecast> triangles;
    Triangulation() = default;
    Triangulation(Delaunator& D);
    int64_t find_triangle_naive(const double x, const double y);
    int64_t find_triangle_bruteforce(const double x, const double y);
    int64_t find_triangle_jump_and_walk(const double x, const double y, uint32_t neighbor);
    bool point_in_triangle(parlay::sequence<double>& bar_coords);
    std::pair<int64_t, int64_t> get_edge(const uint32_t& i, const uint32_t& t);
    std::pair<int64_t, int64_t> check_triangle_line_intersection(uint32_t t, int64_t a, int64_t b, double Ax, double Ay, double Bx, double By);
    std::pair<int64_t, parlay::sequence<double>> jump_and_walk(const double x, const double y, const uint32_t neighbor);
    std::pair<int64_t, parlay::sequence<double>> check_neighborhood(const double x, const double y, const uint32_t pi);
    std::pair<int64_t, parlay::sequence<double>> check_adj_triangles(const double x, const double y, const uint32_t t);
    parlay::sequence<double> barycentric_coordinates(const double x, const double y, const uint32_t t);
    parlay::sequence<double> barycentric_coordinates(const double x, const double y, const uint32_t i, const uint32_t j, const uint32_t k);
};

Triangulation::Triangulation(Delaunator& D) {
//    delauany = &D;
    numVertices = D.coords.size() / 2;
    numTriangles = D.triangles.size() / 3;
    std::vector<int64_t> points_shape = {numVertices, 2};
    std::vector<int64_t> triangles_shape = {numTriangles, 3};
    std::vector<int64_t> points_strides = {sizeof(double) * 2, sizeof(double)};
    std::vector<int64_t> triangles_strides = {sizeof(size_t) * 3, sizeof(size_t)};

    this->vertices = {std::move(points_shape), std::move(points_strides), D.coords.data()};
    this->triangles = {std::move(triangles_shape), std::move(triangles_strides), D.triangles.data()};

    this -> adj_p2t = parlay::sequence<parlay::sequence<uint32_t>>(numVertices);
//    this -> adj_e2t = std::vector<std::vector<std::vector<uint32_t>>>(numVertices, std::vector<std::vector<uint32_t>>(numVertices, std::vector<uint32_t>()));
    uint32_t len = 0;
    for (uint32_t i = 0; i < numTriangles; ++i) {
        for (int32_t j = 0; j < 3; ++j) {
            uint32_t a = triangles.at(i, j);
            uint32_t b = triangles.at(i, fast_mod(j + 1, 3));

            adj_p2t[a].push_back(i);
            if (a > b)
                std::swap(a, b);
            uint64_t e = elegant_pair(a, b);
            auto got = adj_e2t.find(e);
            if (got == adj_e2t.end()) {
                parlay::sequence<uint32_t> s(1, i);
                adj_e2t.insert({e, s});
            }
            else {
                got -> second.push_back(i);
                len = got -> second.size() > len ? got -> second.size() : len;
            }
        }
    }
//    std::cout << edges.size() << ' ' << triangles.shape()[0] << ' ' << vertices.shape()[0] << '\n';
//    std::cout << int32_t(edges.size()) << ' ' << int32_t(triangles.shape()[0]) << ' ' << int32_t(vertices.shape()[0]) << '\n';
//    std::cout << 1 - int32_t(edges.size()) + int32_t(triangles.shape()[0]) + int32_t(vertices.shape()[0]) << '\n';
    std::cout << len << '\n';
}

bool Triangulation::point_in_triangle(parlay::sequence<double>& bar_coords) {
    double r1 = bar_coords[0];
    double r2 = bar_coords[1];
    double r3 = bar_coords[2];
    if ((r1 > 0.0 && r1 < 1.0) and (r2 > 0.0 && r2 < 1.0) and (r3 > 0.0 && r3 < 1.0))
        return true;
    bool mb_on_edge = ((r1 >= 0.0 && r1 <= 1.0) and (r2 >= 0.0 && r2 <= 1.0) and (r3 >= 0.0 && r3 <= 1.0));
    return mb_on_edge and (r1 == 0.0 or r2 == 0.0 or r3 == 0.0);
}

std::pair<int64_t, parlay::sequence<double>> Triangulation::check_adj_triangles(const double x, const double y, const uint32_t t) {
    for (uint32_t i = 0; i < 3; ++i) {
        auto edge = get_edge(i, t);
        uint64_t e = elegant_pair(edge);
        for (auto& tt : adj_e2t[e]) {
            if (tt != t) {
                auto bar_coords = barycentric_coordinates(x, y, tt);
//                std::cout << bar_coords[0].first << ' ' << bar_coords[1].first << ' ' << bar_coords[2].first << std::endl;
                if (point_in_triangle(bar_coords)) {
                    return {tt, bar_coords};
                }
            }
        }
    }
    return {-1, {}};
}

std::pair<int64_t, parlay::sequence<double>> Triangulation::check_neighborhood(const double x, const double y, const uint32_t pi) {
    for (auto& t : this -> adj_p2t[pi]) {   // A bit of overhead :retard:
        auto checked = check_adj_triangles(x, y, t);
        if (checked.first != -1)
            return checked;
    }
    return {-1, {}};
}

parlay::sequence<double> Triangulation::barycentric_coordinates(const double x, const double y, const uint32_t t) {
    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Vertex_approach

    uint32_t i = triangles.at(t, 0);
    uint32_t j = triangles.at(t, 1);
    uint32_t k = triangles.at(t, 2);

    double x1 = vertices.at(i, 0);
    double y1 = vertices.at(i, 1);
    double x2 = vertices.at(j, 0);
    double y2 = vertices.at(j, 1);
    double x3 = vertices.at(k, 0);
    double y3 = vertices.at(k, 1);

    double dx32 = x3 - x2;
    double dx13 = x1 - x3;
    double dx21 = x2 - x1;

    double dy23 = y2 - y3;
    double dy31 = y3 - y1;
    double dy12 = y1 - y2;

    double dxy23 = x2*y3 - x3*y2;
    double dxy31 = x3*y1 - x1*y3;
    double dxy12 = x1*y2 - x2*y1;

    double one_over_2area = 1.0 / double(dxy12 + dxy31 + dxy23);

    auto lambda_1 = double(dxy23 + x * dy23 + y * dx32) * one_over_2area;
    auto lambda_2 = double(dxy31 + x * dy31 + y * dx13) * one_over_2area;
    auto lambda_3 = double(dxy12 + x * dy12 + y * dx21) * one_over_2area;
    return {lambda_1, lambda_2, lambda_3};
}


//parlay::sequence<double> Triangulation::barycentric_coordinates(double x, double y, uint32_t i, uint32_t j, uint32_t k) {
//    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Vertex_approach
//
//    double x1 = vertices.at(i, 0);
//    double y1 = vertices.at(i, 1);
//    double x2 = vertices.at(j, 0);
//    double y2 = vertices.at(j, 1);
//    double x3 = vertices.at(k, 0);
//    double y3 = vertices.at(k, 1);
//
//    double dx32 = x3 - x2;
//    double dx13 = x1 - x3;
//    double dx21 = x2 - x1;
//
//    double dy23 = y2 - y3;
//    double dy31 = y3 - y1;
//    double dy12 = y1 - y2;
//
//    double dxy23 = x2*y3 - x3*y2;
//    double dxy31 = x3*y1 - x1*y3;
//    double dxy12 = x1*y2 - x2*y1;
//
//    double one_over_2area = 1.0 / double(dxy12 + dxy31 + dxy23);
//
//    auto lambda_1 = double(dxy23 + x * dy23 + y * dx32) * one_over_2area;
//    auto lambda_2 = double(dxy31 + x * dy31 + y * dx13) * one_over_2area;
//    auto lambda_3 = double(dxy12 + x * dy12 + y * dx21) * one_over_2area;
//    return {lambda_1, lambda_2, lambda_3};
//}

int64_t Triangulation::find_triangle_naive(const double x, const double y) {
    auto iota = parlay::iota(adj_p2t.size());
    auto minelem = parlay::min_element(iota, [&](auto a, auto b) {
        auto ax = vertices.at(a, 0);
        auto ay = vertices.at(a, 1);
        auto bx = vertices.at(b, 0);
        auto by = vertices.at(b, 1);
        auto dax = (ax - x);
        auto day = (ay - y);
        auto da = dax*dax + day*day;
        auto dbx = (bx - x);
        auto dby = (by - y);
        auto db = dbx*dbx + dby*dby;
        return da < db;
    });
    auto j = minelem - iota.begin();

    return check_neighborhood(x, y, uint32_t(j)).first;
}

int64_t Triangulation::find_triangle_bruteforce(const double x, const double y) {
    for (int64_t i = 0; i < int64_t(numTriangles); ++i) {
        auto bar_coords = barycentric_coordinates(x, y, i);
        bool t = point_in_triangle(bar_coords);
        if (t) {
            return i;
        }
    }
    return -1;
}

std::pair<int64_t, int64_t> Triangulation::get_edge(const uint32_t& i, const uint32_t& t) {
    uint32_t a = triangles.at(t, i);
    uint32_t b = triangles.at(t, fast_mod(i + 1, 3));
    if (a > b)
        return {b, a};
    return {a, b};
}

std::pair<int64_t, int64_t> Triangulation::check_triangle_line_intersection(uint32_t t, int64_t a, int64_t b, double Lx1, double Ly1, double Lx2, double Ly2) {
    for (uint32_t i = 0; i < 3; ++i) {
        auto edge = get_edge(i, t);
        if (edge.first == a and edge.second == b)
            continue;
        double Ax = vertices.at(edge.first, 0), Ay = vertices.at(edge.first, 1);
        double Bx = vertices.at(edge.second, 0), By = vertices.at(edge.second, 1);
        if (segm_intersection(Ax, Ay, Bx, By, Lx1, Ly1, Lx2, Ly2))
            return {edge.first, edge.second};
    }
    return {-1, -1};
}

// Jump-and-Walk algorithm from https://www.cs.montana.edu/bhz/doc/isvd12.pdf

std::pair<int64_t, parlay::sequence<double>> Triangulation::jump_and_walk(const double x, const double y, const uint32_t neighbor) {
    double lx = vertices.at(neighbor, 0);
    double ly = vertices.at(neighbor, 1);
    std::pair<int64_t, int64_t> curr_edge = {-1, -1};
    int64_t curr_tri = -1;
    for (auto& t : adj_p2t[neighbor]) {
        auto bar_coords = barycentric_coordinates(x, y, t);
        if (point_in_triangle(bar_coords)) {
            return {t, bar_coords};
        }
        for (uint32_t i = 0; i < 3; ++i) {
            auto edge = get_edge(i, t);
            if (edge.first == neighbor or edge.second == neighbor)
                continue;
            double ax = vertices.at(edge.first, 0), ay = vertices.at(edge.first, 1);
            double bx = vertices.at(edge.second, 0), by = vertices.at(edge.second, 1);
            if (not segm_intersection(ax, ay, bx, by, x, y, lx, ly))
                continue;
            curr_edge = edge;
            curr_tri = t;
            break;
        }
        if (curr_tri != -1)
            break;
    }

    if (curr_tri == -1)
        return {-1, {}};

    for (uint32_t i = 0; i < 128; ++i) {

        if (curr_edge.first == -1) {
            auto bar_coords = barycentric_coordinates(x, y, curr_tri);
            if (point_in_triangle(bar_coords))
                return {curr_tri, bar_coords};
            return {-1, {}};
        }

        auto adj_tris = adj_e2t[elegant_pair(curr_edge)];
        if (adj_tris[0] == adj_tris[1]) { // check if edge is in dDT
            return {-1, {}};
        }

        //step over current edge

        curr_tri = (adj_tris[0] == curr_tri) ? adj_tris[1] : adj_tris[0];
        auto bar_coords = barycentric_coordinates(x, y, curr_tri);
        if (point_in_triangle(bar_coords))
            return {curr_tri, bar_coords};
        curr_edge = check_triangle_line_intersection(curr_tri, curr_edge.first, curr_edge.second, x, y, lx, ly);

    }

//    int32_t t = find_triangle_bruteforce(x, y);
//    if (t == -1)
//        return {-1, {}};
//    auto bar_coords = barycentric_coordinates(x, y, t);
//
//    return {t, bar_coords};
    curr_tri = find_triangle_bruteforce(x, y);
    if (curr_tri == -1)
        return {-1, {}};
    auto bar_coords = barycentric_coordinates(x, y, curr_tri);
    return {curr_tri, bar_coords};
}

int64_t Triangulation::find_triangle_jump_and_walk(const double x, const double y, const uint32_t neighbor) {
    return jump_and_walk(x, y, neighbor).first;
}

/*
 int32_t tri = -1;
    int32_t e_a = -1;
    int32_t e_b = -1;
    bool find_start = false;
    int32_t lx = vertices.at(neighbor, 0);
    int32_t ly = vertices.at(neighbor, 1);
    for (auto& t : adj_p2t[neighbor]) {
        auto bar_coords = barycentric_coordinates(x, y, t);
        if (point_in_triangle(bar_coords))
            return {t, bar_coords};
        for (uint32_t i = 0; i < 3; ++i) {
            auto edge = get_edge(i, t);
            int32_t a = edge.first;
            int32_t b = edge.second;
            if (a == neighbor or b == neighbor)
                continue;
            int32_t ax = vertices.at(a, 0), ay = vertices.at(a, 1);
            int32_t bx = vertices.at(b, 0), by = vertices.at(b, 1);
            if (segm_intersection(ax, ay, bx, by, x, y, lx, ly)) {
                find_start = true;
//                auto e = bijection(a, b);
                e_a = a;
                e_b = b;
                tri = t;
                break;
            }
        }
        if (find_start)
            break;
    }

    while (true) {
        if (tri == -1)
            return {-1, {}};
        auto bar_coords = barycentric_coordinates(x, y, tri);
        if (point_in_triangle(bar_coords))
            return {tri, bar_coords};
        int32_t prev_tri = tri;
        tri = -1;
        for (auto &t: adj_e2t[bijection(e_a, e_b)]) {
            if (t == prev_tri)
                continue;
            auto inter_edge = check_triangle_line_intersection(t, e_a, e_b, x, y, lx, ly);
            int32_t a = inter_edge.first;
            int32_t b = inter_edge.second;
            int32_t ax = vertices.at(a, 0), ay = vertices.at(a, 1);
            int32_t bx = vertices.at(b, 0), by = vertices.at(b, 1);
            if (segm_intersection(ax, ay, bx, by, x, y, lx, ly)) {
                tri = t;
                e_a = a;
                e_b = b;
                break;
            }
        }
    }
 * */