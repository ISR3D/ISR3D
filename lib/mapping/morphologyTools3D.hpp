#ifndef MORPHOLOGY_TOOLS3D_HPP
#define MORPHOLOGY_TOOLS3D_HPP

#include <vector>

namespace mathematicalMorphologyTools3D {
    bool hasIndirectNeighborWithValue(const std::vector<int32_t>& domain, int i, int j, int k, int32_t value, size_t X, size_t Y, size_t Z) {
        for (int ii = i - 1; ii <= i + 1; ii++) {
            if (ii == -1 || ii == (int) X) continue;
            for (int jj = j - 1; jj <= j + 1; jj++) {
                if (jj == -1 || jj == (int) Y) continue;
                for (int kk = k - 1; kk <= k + 1; kk++) {
                    if (kk == -1 || kk == (int) Z || (ii == i && jj == j && kk == k)) continue;
                    if (domain[ii*Y*Z + jj*Z + kk] == value) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    static void dilate(const std::vector<int32_t>& domain, std::vector<int32_t>& newDomain, size_t X, size_t Y, size_t Z, const int32_t value) {
        for (size_t i = 0; i < X; i++) {
            for (size_t j = 0; j < Y; j++) {
                for (size_t k = 0; k < Z; k++) {
                    if (hasIndirectNeighborWithValue(domain, i, j, k, value, X, Y, Z)) {
                        newDomain[i*Y*Z + j*Z + k] = value;
                    } else {
                        newDomain[i*Y*Z + j*Z + k] = domain[i*Y*Z + j*Z + k];
                    }
                }
            }
        }
    }

    static void erode(const std::vector<int32_t>& domain, std::vector<int32_t>& newDomain, size_t X, size_t Y, size_t Z, const int32_t value) {
        for (size_t i = 0; i < X; i++) {
            for (size_t j = 0; j < Y; j++) {
                for (size_t k = 0; k < Z; k++) {
                    if (domain[i*Y*Z + j*Z + k] == value && hasIndirectNeighborWithValue(domain, i, j, k, 0, X, Y, Z)) {
                        newDomain[i*Y*Z + j*Z + k] = 0;
                    } else {
                        newDomain[i*Y*Z + j*Z + k] = domain[i*Y*Z + j*Z + k];
                    }
                }
            }
        }
    }

    static std::vector<int32_t> open(int times, std::vector<int32_t>& domain, size_t X, size_t Y, size_t Z, const int32_t value) {
        std::vector<int32_t> newDomain(domain.size());

        for (int i = 0; i < times; i++) {
            dilate(domain, newDomain, X, Y, Z, value);
            domain.swap(newDomain);
        }
        for (int i = 0; i < times; i++) {
            erode(domain, newDomain, X, Y, Z, value);
            domain.swap(newDomain);
        }
        return domain;
    }
}

#endif
