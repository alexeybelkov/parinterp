// Szudzik’s Pairing Function https://en.wikipedia.org/wiki/Pairing_function#Other_pairing_functions
inline uint64_t elegant_pair(size_t a, size_t b) {
    return static_cast<uint64_t>(b) * b + a; // implicit cast 
}

inline size_t fast_mod(const size_t i, const size_t c) {
    return i >= c ? i % c : i;
}

// hope it works
inline int64_t ccw(int64_t Ax, int64_t Ay, int64_t Bx, int64_t By, int64_t Cx, int64_t Cy) {
    return (Cy - Ay) * (Bx - Ax) > (By - Ay) * (Cx - Ax);
}

inline bool segments_intersection(size_t Ax, size_t Ay, size_t Bx, size_t By, size_t Cx, size_t Cy, size_t Dx, size_t Dy) {
    return ccw(Ax, Ay, Cx, Cy, Dx, Dy) != ccw(Bx, By, Cx, Cy, Dx, Dy) and ccw(Ax, Ay, Bx, By, Cx, Cy) != ccw(Ax, Ay, Bx, By, Dx, Dy);
}

// hope it also works
inline bool point_in_triangle(std::array<int64_t, 4>& coords_info) {
    return coords_info[3] <= 0 ? 
        coords_info[3] <= coords_info[0] and coords_info[0] <= 0
        and coords_info[3] <= coords_info[1] and coords_info[1] <= 0 
        and coords_info[3] <= coords_info[2] and coords_info[2] <= 0 : 
        coords_info[3] >= coords_info[0] and coords_info[0] >= 0
        and coords_info[3] >= coords_info[1] and coords_info[1] >= 0 
        and coords_info[3] >= coords_info[2] and coords_info[2] >= 0;
}
