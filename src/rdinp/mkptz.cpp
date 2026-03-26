// Build polarization tensor and rotate coordinates.
// Converted from: src/RDINP/mkptz.f

#include "mkptz.hpp"
#include <feff/constants.hpp>
#include "../common/logging.hpp"

#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace feff::rdinp {

void rotate(double vec[3], double cst, double snt, double csf, double snf) {
    double temp[3];
    temp[0] = vec[0] * cst * csf + vec[1] * cst * snf - vec[2] * snt;
    temp[1] = -vec[0] * snf + vec[1] * csf;
    temp[2] = vec[0] * csf * snt + vec[1] * snt * snf + vec[2] * cst;
    vec[0] = temp[0];
    vec[1] = temp[1];
    vec[2] = temp[2];
}

void mkptz(int ipol, double elpty, double evec[3], double xivec[3],
           int ispin, double spvec[3], int nat, double rat[][3],
           double& angks, int& le2, FeffComplex ptz[3][3]) {

    auto& log = feff::common::logger();

    // Make z axis along propagation (XIVEC).
    double rr = xivec[0] * xivec[0] + xivec[1] * xivec[1] + xivec[2] * xivec[2];
    if (rr == 0.0) {
        angks = 0.0;
        // Special case: xivec not specified
        if (ipol == 1) {
            // Need xivec for E2 and M1 transitions; leave only E1
            if (le2 != 0) {
                log.wlog("  Can do only E1 transitions. Specify k-vector for M1 or E2");
            }
            le2 = 0;
        } else {
            // For polarization average or circular dichroism
            if (ispin != 0) {
                // Spin-dependent case: use spvec as xivec
                for (int i = 0; i < 3; ++i) {
                    xivec[i] = spvec[i];
                }
                rr = xivec[0] * xivec[0] + xivec[1] * xivec[1] + xivec[2] * xivec[2];
            }
        }
    }

    if (rr > 0.0) {
        double rsp = std::sqrt(rr);
        rr = xivec[0] * xivec[0] + xivec[1] * xivec[1];
        if (rr != 0.0 || xivec[2] < 0.0) {
            double cst, snt, csf, snf;
            if (rr == 0.0) {
                cst = -1.0;
                snt = 0.0;
                csf = 1.0;
                snf = 0.0;
            } else {
                rr = std::sqrt(rr);
                cst = xivec[2] / rsp;
                snt = rr / rsp;
                csf = xivec[0] / rr;
                snf = xivec[1] / rr;
            }
            // Rotate all vectors
            for (int i = 0; i < nat; ++i) {
                rotate(rat[i], cst, snt, csf, snf);
            }
            rotate(evec, cst, snt, csf, snf);
            rotate(xivec, cst, snt, csf, snf);
            rotate(spvec, cst, snt, csf, snf);
        }
    }

    // Initialize ptz to zero
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ptz[i][j] = FeffComplex(0.0, 0.0);
        }
    }

    // Build ptz in frame where z is along xivec (except ipol=0)
    if (ipol == 0) {
        // Orientation average: ptz(i,i) = 1/3
        for (int i = 0; i < 3; ++i) {
            ptz[i][i] = FeffComplex(1.0 / 3.0, 0.0);
        }
    } else if (ipol == 2) {
        // Circular polarization
        ptz[2][2] = FeffComplex(1.0, 0.0);    // ptz(1,1) in Fortran (-1:1 indexing)
        ptz[0][0] = FeffComplex(-1.0, 0.0);   // ptz(-1,-1) in Fortran
    } else if (ipol == 1) {
        // Linear/elliptical polarization

        // Normalize polarization vector
        double x = std::sqrt(evec[0] * evec[0] + evec[1] * evec[1] + evec[2] * evec[2]);
        if (x <= 0.000001) {
            log.wlog(" STOP  Polarization vector of almost zero length");
            log.wlog(" Correct POLARIZATION card");
            throw std::runtime_error("MKPTZ-1: zero-length polarization vector");
        }
        for (int i = 0; i < 3; ++i) {
            evec[i] /= x;
        }

        x = std::sqrt(xivec[0] * xivec[0] + xivec[1] * xivec[1] + xivec[2] * xivec[2]);
        if (x > 0.0) {
            // Run elliptical polarization code
            for (int i = 0; i < 3; ++i) {
                xivec[i] /= x;
            }

            // Check that evec is not parallel to xivec
            x = evec[0] * xivec[0] + evec[1] * xivec[1] + evec[2] * xivec[2];
            if (std::abs(x) > 0.9) {
                std::ostringstream ss;
                log.wlog(" polarization");
                ss << "     " << std::scientific << std::setprecision(5)
                   << evec[0] << "  " << evec[1] << "  " << evec[2];
                log.wlog(ss.str());
                ss.str(""); ss.clear();
                log.wlog(" incidence");
                ss << "     " << std::scientific << std::setprecision(5)
                   << xivec[0] << "  " << xivec[1] << "  " << xivec[2];
                log.wlog(ss.str());
                log.wlog(" STOP polarization almost parallel to the incidence");
                log.wlog(" Correct ELLIPTICITY and POLARIZATION cards");
                throw std::runtime_error("MKPTZ-2: polarization parallel to incidence");
            }

            if (x != 0.0) {
                // Make evec normal to xivec, keeping the plane
                log.wlog(" Changing polarization vector!");
                log.wlog(" Incidence is not normal to polarization.");
                log.wlog(" Check your input for errors. Run continues.");
                for (int i = 0; i < 3; ++i) {
                    evec[i] = evec[i] - x * xivec[i];
                }
                x = std::sqrt(evec[0] * evec[0] + evec[1] * evec[1] + evec[2] * evec[2]);
                for (int i = 0; i < 3; ++i) {
                    evec[i] /= x;
                }
            }
        } else {
            // elpty cannot be used with xivec=0
            elpty = 0.0;
        }

        // e2 = xivec cross evec
        double e2[3];
        e2[0] = xivec[1] * evec[2] - xivec[2] * evec[1];
        e2[1] = xivec[2] * evec[0] - xivec[0] * evec[2];
        e2[2] = xivec[0] * evec[1] - xivec[1] * evec[0];

        // Build complex polarization vector e = evec + i*elpty*e2
        FeffComplex e[3], eps_arr[3], epc_arr[3]; // eps(-1:1), epc(-1:1) -> [0]=(-1),[1]=(0),[2]=(1)
        for (int i = 0; i < 3; ++i) {
            e[i] = FeffComplex(evec[i], 0.0) + elpty * FeffComplex(0.0, e2[i]);
        }
        // eps(-1) = (e(1) - i*e(2)) / sqrt(2)
        eps_arr[0] = (e[0] - coni * e[1]) / std::sqrt(2.0);
        // eps(0) = e(3)
        eps_arr[1] = e[2];
        // eps(1) = -(e(1) + i*e(2)) / sqrt(2)
        eps_arr[2] = -(e[0] + coni * e[1]) / std::sqrt(2.0);

        // Build complex conjugate polarization: epc = evec - i*elpty*e2
        for (int i = 0; i < 3; ++i) {
            e[i] = FeffComplex(evec[i], 0.0) - elpty * FeffComplex(0.0, e2[i]);
        }
        epc_arr[0] = (e[0] - coni * e[1]) / std::sqrt(2.0);
        epc_arr[1] = e[2];
        epc_arr[2] = -(e[0] + coni * e[1]) / std::sqrt(2.0);

        // Build polarization tensor
        // Fortran: ptz(j,i) for i,j = -1,0,1
        // C++: ptz[j+1][i+1] to match Fortran's ptz(j,i) storage order
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                double sign = (i % 2 == 0) ? 1.0 : -1.0; // (-1)^i
                // Average over left and right for path reversal symmetry
                ptz[j + 1][i + 1] = sign *
                    (epc_arr[i + 1] * eps_arr[-j + 1] + eps_arr[i + 1] * epc_arr[-j + 1])
                    / (1.0 + elpty * elpty) / 2.0;
            }
        }
    }
    // End of making polarization tensor

    angks = 0.0;

    // Second rotation: make z parallel to spin direction
    rr = spvec[0] * spvec[0] + spvec[1] * spvec[1] + spvec[2] * spvec[2];
    if (rr > 0.0) {
        double rsp = std::sqrt(rr);
        rr = spvec[0] * spvec[0] + spvec[1] * spvec[1];
        if (rr != 0.0 || spvec[2] < 0.0) {
            double cst, snt, csf, snf;
            if (rr == 0.0) {
                cst = -1.0;
                snt = 0.0;
                csf = 1.0;
                snf = 0.0;
                angks = pi;
            } else {
                rr = std::sqrt(rr);
                cst = spvec[2] / rsp;
                snt = rr / rsp;
                csf = spvec[0] / rr;
                snf = spvec[1] / rr;
                angks = std::acos(cst);
            }
            // Rotate all vectors
            for (int i = 0; i < nat; ++i) {
                rotate(rat[i], cst, snt, csf, snf);
            }
            rotate(evec, cst, snt, csf, snf);
            rotate(xivec, cst, snt, csf, snf);
        }
    }
}

} // namespace feff::rdinp
