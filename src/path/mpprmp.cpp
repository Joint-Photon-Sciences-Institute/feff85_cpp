// Make path parameters xp, yp, zp in standard reference frame.
// Converted from: src/PATH/mpprmp.f

#include "mpprmp.hpp"
#include <cmath>

namespace feff::path {

void mpprmp(int npat, const int ipat[], float xp[], float yp[], float zp[],
            int ipol, int ispin, const double evec[3], const double xivec[3],
            int ica, const AtomData& atoms) {

    constexpr float eps4 = 1.0e-4f;

    int nleg = npat + 1;

    // Determine symmetry case
    int icase;
    if (ica > 0 && ica < 8) {
        icase = ica;
    } else {
        bool lkvec = (xivec[0]*xivec[0] + xivec[1]*xivec[1] + xivec[2]*xivec[2]) > eps4;
        bool lkz   = (xivec[0]*xivec[0] + xivec[1]*xivec[1]) <= eps4;
        bool lez   = (evec[0]*evec[0] + evec[1]*evec[1]) <= eps4;

        icase = 7; // default: no symmetry
        if (ipol == 0) {
            icase = 1;
        } else if (ispin == 0) {
            if (ipol == 1 && !lkvec)  icase = 2;
            if (ipol == 1 && lkvec)   icase = 3;
            if (ipol == 2)            icase = 4;
        } else {
            if (ipol == 2 && lkz)     icase = 5;
            if (ipol == 1 && !lkvec && lez) icase = 5;
            if (ipol == 1 && lkz && evec[2]*evec[2] < eps4) icase = 6;
        }
    }

    // Initialize outputs
    std::vector<double> xp1(npatx, 0.0), yp1(npatx, 0.0), zp1(npatx, 0.0);

    for (int j = 0; j < npatx; ++j) {
        xp[j] = 0.0f;
        yp[j] = 0.0f;
        zp[j] = 0.0f;
    }

    // ri[i][j] = rat(i, ipat[j]) - rat(i, 0) for j=0..npat-1
    double ri[3][npatx];
    for (int j = 0; j < npat; ++j) {
        for (int i = 0; i < 3; ++i) {
            ri[i][j] = atoms.rat[ipat[j]][i] - atoms.rat[0][i];
        }
    }
    for (int j = nleg - 1; j < npatx; ++j) {
        for (int i = 0; i < 3; ++i) {
            ri[i][j] = 0.0;
        }
    }

    double xvec[3] = {0, 0, 0};
    double yvec[3] = {0, 0, 0};
    double zvec[3] = {0, 0, 0};

    // Choose z-axis
    if (icase == 1) {
        double norm = 0.0;
        for (int i = 0; i < 3; ++i) norm += ri[i][0] * ri[i][0];
        norm = std::sqrt(norm);
        for (int i = 0; i < 3; ++i) zvec[i] = ri[i][0] / norm;
    } else if (icase == 2 || icase == 3) {
        for (int i = 0; i < 3; ++i) zvec[i] = evec[i];
    } else {
        zvec[2] = 1.0;
    }

    // Compute z-coordinates
    for (int j = 0; j < npat; ++j) {
        zp1[j] = 0.0;
        for (int i = 0; i < 3; ++i) {
            zp1[j] += zvec[i] * ri[i][j];
        }
    }

    // If no symmetries, don't waste time
    if (icase == 7) {
        xvec[0] = 1.0;
        yvec[1] = 1.0;
        goto label_390;
    }

    // Use z --> -z symmetry (cases 2, 3, 6)
    if (!(icase == 1 || icase >= 4)) {
        for (int num = 0; num < nleg - 1; ++num) {
            if (std::abs(zp1[num]) > eps4) {
                if (zp1[num] < 0) {
                    for (int j = 0; j < 3; ++j) zvec[j] = -zvec[j];
                    for (int j = 0; j < npat; ++j) zp1[j] = -zp1[j];
                }
                goto label_240;
            }
        }
    }

label_240:
    // Use rotations around z and reflections containing z
    for (int num = 0; num < nleg - 1; ++num) {
        double ro2 = 0.0;
        for (int i = 0; i < 3; ++i) ro2 += ri[i][num] * ri[i][num];
        ro2 = ro2 - zp1[num] * zp1[num];
        ro2 = std::sqrt(std::abs(ro2));

        if (ro2 >= eps4) {
            // Atom not on z-axis
            if (icase == 1 || icase == 2 || icase == 4 || icase == 5) {
                // Any rotation around z: choose x so that x-coord positive, y=0
                for (int i = 0; i < 3; ++i) {
                    xvec[i] = ri[i][num] - zvec[i] * zp1[num];
                }
                for (int i = 0; i < 3; ++i) xvec[i] /= ro2;
            } else if (icase == 3) {
                // Elliptical polarization: x along incident beam
                for (int i = 0; i < 3; ++i) xvec[i] = xivec[i];
            } else {
                // icase == 6: choose x so that x-coord is positive
                xvec[0] = 1.0;
                if (ri[0][num] < 0) xvec[0] = -1.0;
            }
            // y = z cross x
            yvec[0] = zvec[1] * xvec[2] - zvec[2] * xvec[1];
            yvec[1] = zvec[2] * xvec[0] - zvec[0] * xvec[2];
            yvec[2] = zvec[0] * xvec[1] - zvec[1] * xvec[0];
            goto label_390;
        }
    }

label_390:
    // Calculate x,y coords in the chosen frame
    for (int j = 0; j < npat; ++j) {
        xp1[j] = 0.0;
        yp1[j] = 0.0;
        for (int i = 0; i < 3; ++i) {
            xp1[j] += xvec[i] * ri[i][j];
            yp1[j] += yvec[i] * ri[i][j];
        }
    }

    // For icase==3: check first nonzero x is positive
    if (icase == 3) {
        for (int num = 0; num < nleg - 1; ++num) {
            if (std::abs(xp1[num]) >= eps4) {
                if (xp1[num] < 0) {
                    for (int j = 0; j < npat; ++j) xp1[j] = -xp1[j];
                }
                break;
            }
        }
    }

    // For cases 1-3: inverse y if first nonzero y is negative
    if (icase < 4) {
        for (int num = 0; num < nleg - 1; ++num) {
            if (std::abs(yp1[num]) >= eps4) {
                if (yp1[num] < 0) {
                    for (int j = 0; j < npat; ++j) yp1[j] = -yp1[j];
                }
                break;
            }
        }
    }

    // Copy to output (single precision)
    for (int j = 0; j < npat; ++j) {
        xp[j] = static_cast<float>(xp1[j]);
        yp[j] = static_cast<float>(yp1[j]);
        zp[j] = static_cast<float>(zp1[j]);
    }
}

} // namespace feff::path
