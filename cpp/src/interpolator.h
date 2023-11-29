// using namespace py::literals;
class Interpolator {
private:
    size_t n_jobs_;
public:
    using pyarr_double = py::array_t<double, py::array::c_style | py::array::forcecast>;
    const Triangulator triangulation;
    py::object kdtree;
    Interpolator(const Triangulator::pyarr_size_t& points, size_t n_jobs) : triangulation(points, n_jobs), n_jobs_(n_jobs) {
        // kdtree = py::module_::import("scipy").attr("spatial").attr("KDTree")(points);
    }
    pyarr_double operator()(const Triangulator::pyarr_size_t& int_points, const pyarr_double& values, const Triangulator::pyarr_size_t& neighbors, double fill_value = 0.0) {
        // using namespace py::literals;
        if (triangulation.ppoints -> shape()[0] != values.shape()[0]) {
            throw std::invalid_argument("Length mismath between known points and their values");
        }
        if (neighbors.shape()[0] != int_points.shape()[0]) {
            throw std::invalid_argument("Length mismath between int_points and their neighbors");
        }
        size_t n = int_points.shape()[0];
        std::vector<double> int_values(n, fill_value);
        // auto neighbors = kdtree.attr("query")(int_points, "workers"_a = n_jobs_)["1"];
        // omp_set_num_threads(n_jobs_);
        // #pragma parallel for
        for (size_t i = 0; i < n; ++i) {
            size_t x = int_points.at(i, 0);
            size_t y = int_points.at(i, 1);
            auto location_info = triangulation.locate_point(x, y, neighbors.at(i));
            if (location_info != std::nullopt) {
                size_t t = location_info -> first;
                auto& coords_info = location_info -> second;
                double one_over_2area = 1.0 / static_cast<double>(coords_info[3]);
                double lambda_1 = static_cast<double>(coords_info[0]) / one_over_2area;
                double lambda_2 = static_cast<double>(coords_info[1]) / one_over_2area;
                double lambda_3 = static_cast<double>(coords_info[2]) / one_over_2area;
                double f1 = values.at(triangulation.triangles.at(t));
                double f2 = values.at(triangulation.triangles.at(t + 1));
                double f3 = values.at(triangulation.triangles.at(t + 2));
                int_values[i] = f1 * lambda_1 + f2 * lambda_2 + f3 * lambda_3;
            }
        }
        return {{n}, {sizeof(double)}, int_values.data()};
    }
};
