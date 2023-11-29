//
// Created by lbelkov on 8/31/23.
//

#ifndef PARDELAUNAY_RADIXSORT_H
#define PARDELAUNAY_RADIXSORT_H

#endif //PARDELAUNAY_RADIXSORT_H

#include <vector>
#include <cstdint>

void radix_sort(std::vector<std::size_t>& data, const size_t BITS = 8) {
    const size_t BUCKET = 1 << BITS;
    const size_t LAYERS = sizeof(size_t) / BITS;

    std:: vector<size_t> cnt(BUCKET);
    std:: vector<size_t> str(BUCKET);
    std:: vector<size_t> end(BUCKET);
    std::vector<std::size_t> copy(data.size());

    for (int layer = 0; layer < int(LAYERS); ++layer) {
        cnt.assign(BUCKET, 0);
        for (auto it = data.begin(); it != data.end(); ++it) {
            std::uint64_t key = (*it >> (BITS * layer)) & (BUCKET - 1);
            ++cnt[key];
        }
        str[0] = 0;
        end[0] = cnt[0];
        for (size_t i = 1; i < BUCKET; ++i) {
            str[i] = str[i - 1] + cnt[i - 1];
            end[i] = end[i - 1] + cnt[i];
        }
        for (auto it = data.begin(); it != data.end(); ++it) {
            uint64_t key = (*it >> (BITS * layer)) & (BUCKET - 1);
            copy[str[key]++] = *it;
        }
        data.swap(copy);
    }
}
