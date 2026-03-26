// BiCGStab solver for FMS.
// Converted from: ggbi.f
// Inverts (1 - T*G0) column-by-column using BiCGStab iteration.

#include "ggbi.hpp"

#include <cmath>
#include <complex>
#include <iostream>

namespace feff::fms {

// Helper: complex dot product <bra|ket> with conjugation of bra
static Complexf cdot(int n, const Complexf* bra, const Complexf* ket) {
    Complexf sum(0.0f, 0.0f);
    for (int i = 0; i < n; ++i) {
        sum += std::conj(bra[i]) * ket[i];
    }
    return sum;
}

// Helper: check if all elements have |re| and |im| < tol
static bool converged(int n, const Complexf* v, float tol) {
    for (int i = 0; i < n; ++i) {
        if (std::abs(v[i].real()) > tol || std::abs(v[i].imag()) > tol)
            return false;
    }
    return true;
}

void ggbi(int nsp, const int* i0, int ipi, int ipf, const int* lipotx,
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

    // Build g0t = I - T*G0  (note: ggbi uses T*G0, not G0*T like gglu)
    g0t = Eigen::MatrixXcf::Zero(istate, istate);

    for (int icol = 0; icol < istate; ++icol) {
        for (int irow = 0; irow < istate; ++irow) {
            if (std::abs(g0(irow, icol)) > toler2) {
                g0t(irow, icol) -= tmatrx_diag(irow) * g0(irow, icol);
            }
            // Spin-flip contribution
            int l1   = basis.l(irow);
            int m1   = basis.m(irow);
            int isp1 = basis.isp(irow);
            int m2   = m1 + isp1;
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

    // Solve column-by-column for each potential and l-channel
    Eigen::VectorXcf xvec(istate), rvec(istate), pvec(istate);
    Eigen::VectorXcf svec(istate), avec(istate), asve(istate), yvec(istate);

    for (int ip = ipi; ip <= ipf; ++ip) {
        int ipart = nsp * (lipotx[ip] + 1) * (lipotx[ip] + 1);
        for (int is1 = 0; is1 < ipart; ++is1) {
            int is2 = is1 + i0[ip];
            int l1 = basis.l(is2);
            if (!lcalc[l1]) continue;

            // BiCGStab iteration to solve g0t * xvec = e_{is2}
            xvec.setZero();
            int istart = -1;
            int local_msord = 0;

        restart:
            ++istart;

            if (istart > 0) {
                avec.noalias() = g0t.topLeftCorner(istate, istate) * xvec;
                local_msord++;
            } else {
                avec.setZero();
            }

            for (int is = 0; is < istate; ++is)
                rvec(is) = -avec(is);
            rvec(is2) += Complexf(1.0f, 0.0f);

            if (converged(istate, rvec.data(), toler1)) goto done;

            pvec = rvec;
            avec.noalias() = g0t.topLeftCorner(istate, istate) * pvec;
            local_msord++;

            // Choose yvec for good conditioning
            {
                Complexf aa = cdot(istate, avec.data(), avec.data());
                Complexf wa = cdot(istate, rvec.data(), avec.data());
                Complexf aw = std::conj(wa);
                Complexf ww = cdot(istate, rvec.data(), rvec.data());
                Complexf dd = aa * ww - aw * wa;
                if (std::abs(dd / (aa * ww)) < 1.0e-8f) {
                    yvec = rvec / ww;
                } else {
                    Complexf ww2 = (ww - aw) / dd;
                    Complexf aa2 = (wa - aa) / dd;
                    yvec = rvec * aa2 + avec * ww2;
                }
            }

            {
                Complexf del = cdot(istate, yvec.data(), rvec.data());
                constexpr int nitx = 30;

                for (int nit = 0; nit <= nitx; ++nit) {
                    Complexf delp = cdot(istate, yvec.data(), avec.data());
                    Complexf omega = del / delp;

                    svec = rvec - omega * avec;
                    if (converged(istate, svec.data(), toler1)) {
                        xvec += omega * pvec;
                        goto done;
                    }

                    asve.noalias() = g0t.topLeftCorner(istate, istate) * svec;
                    local_msord++;

                    Complexf aa2 = cdot(istate, asve.data(), asve.data());
                    Complexf wa2 = cdot(istate, asve.data(), svec.data());
                    Complexf chi = wa2 / aa2;

                    xvec += omega * pvec + chi * svec;
                    rvec = svec - chi * asve;

                    if (converged(istate, rvec.data(), toler1)) goto done;

                    // Prepare next iteration
                    del = cdot(istate, yvec.data(), rvec.data());
                    Complexf psi = del / (delp * chi);
                    pvec = rvec + psi * (pvec - chi * avec);
                    avec.noalias() = g0t.topLeftCorner(istate, istate) * pvec;
                    local_msord++;
                }
                // Ran out of iterations, restart
                goto restart;
            }

        done:
            msord += local_msord;

            // Pack result: gg(is2_out, is1, ip) = sum_is g0(is2_out+i0[ip], is) * xvec(is)
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
