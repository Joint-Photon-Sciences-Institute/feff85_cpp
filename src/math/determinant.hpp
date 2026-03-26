#pragma once

// Matrix determinant calculation
// Converted from src/MATH/determ.f

namespace feff::math {

// Calculate determinant of a square matrix using Gaussian elimination.
// array: matrix of size nrows x nrows (row-major). Modified in place.
// nord: order of matrix (active submatrix size).
// nrows: leading dimension of array.
double determ(double* array, int nord, int nrows);

} // namespace feff::math
