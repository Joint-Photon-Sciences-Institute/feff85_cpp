#include "determinant.hpp"
#include <cmath>

namespace feff::math {

// Access element (i,j) of a matrix stored in row-major with leading dim nrows
static inline double& mat(double* array, int i, int j, int nrows) {
    return array[i * nrows + j];
}

double determ(double* array, int nord, int nrows) {
    double det = 1.0;

    for (int k = 0; k < nord; ++k) {
        if (mat(array, k, k, nrows) == 0.0) {
            // Find non-zero pivot in row k
            int j;
            for (j = k; j < nord; ++j) {
                if (mat(array, k, j, nrows) != 0.0) break;
            }
            if (j >= nord) {
                return 0.0;  // Singular
            }
            // Swap columns k and j
            for (int i = k; i < nord; ++i) {
                double saved = mat(array, i, j, nrows);
                mat(array, i, j, nrows) = mat(array, i, k, nrows);
                mat(array, i, k, nrows) = saved;
            }
            det = -det;
        }

        det *= mat(array, k, k, nrows);

        if (k < nord - 1) {
            for (int i = k + 1; i < nord; ++i) {
                for (int jj = k + 1; jj < nord; ++jj) {
                    mat(array, i, jj, nrows) -= mat(array, i, k, nrows) *
                                                 mat(array, k, jj, nrows) /
                                                 mat(array, k, k, nrows);
                }
            }
        }
    }

    return det;
}

} // namespace feff::math
