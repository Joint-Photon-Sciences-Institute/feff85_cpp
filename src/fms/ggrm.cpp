// Recursion method (non-symmetric Lanczos) solver for FMS.
// Converted from: ggrm.f
// Includes helper functions: cdot, matvec (also used by ggbi, gggm, ggtf).

#include "ggrm.hpp"

#include <cmath>
#include <complex>

namespace feff::fms {

// Helper: conjugate dot product <bra*|ket>
static Complexf cdot_rm(int n, const Complexf* bra, const Complexf* ket) {
    Complexf sum(0.0f, 0.0f);
    for (int i = 0; i < n; ++i) {
        sum += std::conj(bra[i]) * ket[i];
    }
    return sum;
}

static bool converged_rm(int n, const Complexf* v, float tol) {
    for (int i = 0; i < n; ++i) {
        if (std::abs(v[i].real()) > tol || std::abs(v[i].imag()) > tol)
            return false;
    }
    return true;
}

// Build g0t = I - T*G0 (same pattern as ggbi)
static void build_g0t_tg(int nsp, int istate,
                          const Eigen::MatrixXcf& g0,
                          const Eigen::VectorXcf& tmatrx_diag,
                          const Eigen::VectorXcf& tmatrx_off,
                          Eigen::MatrixXcf& g0t,
                          float toler2,
                          const BasisStates& basis) {
    g0t = Eigen::MatrixXcf::Zero(istate, istate);
    for (int icol = 0; icol < istate; ++icol) {
        for (int irow = 0; irow < istate; ++irow) {
            if (std::abs(g0(irow, icol)) > toler2) {
                g0t(irow, icol) -= tmatrx_diag(irow) * g0(irow, icol);
            }
            int l1 = basis.l(irow);
            int m1 = basis.m(irow);
            int isp1 = basis.isp(irow);
            int m2 = m1 + isp1;
            if (nsp == 2 && m2 > -l1 + 1 && m2 < l1 + 2) {
                int ist2 = irow + ((isp1 == 0) ? 1 : -1);
                if (ist2 >= 0 && ist2 < istate &&
                    std::abs(g0(ist2, icol)) > toler2) {
                    g0t(irow, icol) -= tmatrx_off(ist2) * g0(ist2, icol);
                }
            }
        }
        g0t(icol, icol) += Complexf(1.0f, 0.0f);
    }
}

void ggrm(int nsp, const int* i0, int ipi, int ipf, const int* lipotx,
          Eigen::MatrixXcf& g0,
          const Eigen::VectorXcf& tmatrx_diag,
          const Eigen::VectorXcf& tmatrx_off,
          Eigen::MatrixXcf& g0t,
          std::vector<Eigen::MatrixXcf>& gg,
          float toler1, float toler2, const bool* lcalc, int& msord,
          const FMSData& data) {

    const auto& basis = data.basis;
    int istate = basis.istate;
    msord = 0;

    build_g0t_tg(nsp, istate, g0, tmatrx_diag, tmatrx_off, g0t, toler2, basis);

    Eigen::VectorXcf xvec(istate), xket(istate), xbra(istate);
    Eigen::VectorXcf xketp(istate), xbrap(istate);
    Eigen::VectorXcf zvec(istate), rvec(istate), svec(istate);
    Eigen::VectorXcf tket(istate), tbra(istate), avec(istate);

    Complexf coni_f(0.0f, 1.0f);

    for (int ip = ipi; ip <= ipf; ++ip) {
        int ipart = nsp * (lipotx[ip] + 1) * (lipotx[ip] + 1);
        for (int is1 = 0; is1 < ipart; ++is1) {
            int is2 = is1 + i0[ip];
            int l1 = basis.l(is2);
            if (!lcalc[l1]) continue;

            xvec.setZero();
            rvec.setZero();
            int istart = -1;
            int local_msord = 0;

        restart:
            ++istart;

            if (istart > 0) {
                rvec.noalias() = g0t.topLeftCorner(istate, istate) * xvec;
                local_msord++;
            }
            // rvec = g0t*xvec - bvec
            rvec(is2) -= Complexf(1.0f, 0.0f);
            xket = -rvec;

            {
                Complexf bb = cdot_rm(istate, xket.data(), xket.data());
                if (std::abs(bb) == 0.0f) goto done;

                float xfnorm = 1.0f / bb.real();
                xbra = xket * xfnorm;

                // |t> = A|n>
                tket.noalias() = g0t.topLeftCorner(istate, istate) * xket;
                local_msord++;

                Complexf aa = cdot_rm(istate, xbra.data(), tket.data());
                Complexf aac = std::conj(aa);
                bb = Complexf(0.0f, 0.0f);
                Complexf bbc(0.0f, 0.0f);
                Complexf betac = aa;
                Complexf yy(1.0f, 0.0f);

                xketp.setZero();
                xbrap.setZero();
                zvec = xket;
                xvec += zvec / betac;
                svec = tket;
                rvec += svec / betac;

                constexpr int nitx = 100;
                for (int nit = 1; nit <= nitx; ++nit) {
                    // Recursion: update n+1 Lanczos vectors
                    tket = tket - aa * xket - bb * xketp;
                    // Adjoint recursion
                    tbra.noalias() = g0t.topLeftCorner(istate, istate).adjoint() * xbra;
                    tbra = tbra - aac * xbra - bbc * xbrap;

                    bb = cdot_rm(istate, tbra.data(), tket.data());
                    if (std::abs(bb) == 0.0f) goto done;

                    bb = std::sqrt(bb);
                    bbc = std::conj(bb);

                    xketp = xket;
                    xbrap = xbra;
                    xket = tket / bb;
                    xbra = tbra / bbc;

                    tket.noalias() = g0t.topLeftCorner(istate, istate) * xket;
                    local_msord++;

                    aa = cdot_rm(istate, xbra.data(), tket.data());
                    aac = std::conj(aa);

                    // Update solution
                    Complexf alphac = bb / betac;
                    zvec = xket - alphac * zvec;
                    svec = tket - alphac * svec;
                    betac = aa - alphac * bb;
                    yy = -alphac * yy;
                    Complexf gamma = yy / betac;

                    xvec += gamma * zvec;
                    rvec += gamma * svec;

                    if (converged_rm(istate, rvec.data(), toler1)) goto done;
                }
                // Restart
                goto restart;
            }

        done:
            msord += local_msord;

            // Pack into gg
            for (int is2_out = 0; is2_out < ipart; ++is2_out) {
                Complexf val(0.0f, 0.0f);
                for (int is = 0; is < istate; ++is) {
                    val += g0(is2_out + i0[ip], is) * xvec(is);
                }
                gg[ip](is2_out, is1) = val;
            }
        }
    }
}

} // namespace feff::fms
