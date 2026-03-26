// TFQMR solver for FMS.
// Converted from: ggtf.f
// Reference: Saad, "Iterative Methods for Sparse Matrices", p.225 (1996)

#include "ggtf.hpp"

#include <cmath>
#include <complex>

namespace feff::fms {

static Complexf cdot_tf(int n, const Complexf* bra, const Complexf* ket) {
    Complexf sum(0.0f, 0.0f);
    for (int i = 0; i < n; ++i)
        sum += std::conj(bra[i]) * ket[i];
    return sum;
}

void ggtf(int nsp, const int* i0, int ipi, int ipf, const int* lipotx,
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

    // Build g0t = I - T*G0
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

    Eigen::VectorXcf xvec(istate), uvec(istate), avec(istate), wvec(istate);
    Eigen::VectorXcf dvec(istate), rvec(istate), vvec(istate);

    for (int ip = ipi; ip <= ipf; ++ip) {
        int ipart = nsp * (lipotx[ip] + 1) * (lipotx[ip] + 1);
        for (int is1 = 0; is1 < ipart; ++is1) {
            int is2 = is1 + i0[ip];
            int l1_check = basis.l(is2);
            if (!lcalc[l1_check]) continue;

            xvec.setZero();
            rvec.setZero();
            avec.setZero();
            int local_msord = 0;
            bool solved = false;

            Complexf alpha(1.0f, 0.0f);

            for (int istart = 0; !solved; ++istart) {
                if (istart > 0) {
                    avec.noalias() = g0t.topLeftCorner(istate, istate) * xvec;
                    local_msord++;
                }

                // uvec = bvec - g0t*xvec
                uvec = -avec;
                uvec(is2) += Complexf(1.0f, 0.0f);

                avec.noalias() = g0t.topLeftCorner(istate, istate) * uvec;
                local_msord++;

                wvec = uvec;
                vvec = avec;
                dvec.setZero();

                Complexf aa = cdot_tf(istate, uvec.data(), uvec.data());
                float tau = std::sqrt(aa.real());
                float nu = 0.0f;
                Complexf eta(0.0f, 0.0f);

                // Choose rvec = uvec / aa
                rvec = uvec / aa;
                Complexf rho(1.0f, 0.0f);

                constexpr int nitx = 20;
                for (int nit = 0; nit <= nitx; ++nit) {
                    if (nit % 2 == 0) {
                        aa = cdot_tf(istate, rvec.data(), vvec.data());
                        alpha = rho / aa;
                    } else {
                        avec.noalias() = g0t.topLeftCorner(istate, istate) * uvec;
                        local_msord++;
                    }

                    wvec -= alpha * avec;

                    aa = Complexf(nu * nu, 0.0f) * eta / alpha;
                    dvec = uvec + aa * dvec;

                    Complexf ww = cdot_tf(istate, wvec.data(), wvec.data());
                    nu = std::sqrt(ww.real()) / tau;
                    float cm = 1.0f / std::sqrt(1.0f + nu * nu);
                    tau = tau * nu * cm;
                    eta = Complexf(cm * cm, 0.0f) * alpha;

                    xvec += eta * dvec;

                    // Convergence check
                    float err = static_cast<float>(1 + nit) / istate;
                    err = tau * std::sqrt(err) * 10.0f;
                    if (std::abs(err) < toler1) {
                        solved = true;
                        break;
                    }

                    if (nit % 2 != 0) {
                        aa = rho;
                        rho = cdot_tf(istate, rvec.data(), wvec.data());
                        Complexf beta = rho / aa;
                        uvec = wvec + beta * uvec;
                        vvec = beta * (avec + beta * vvec);
                        avec.noalias() = g0t.topLeftCorner(istate, istate) * uvec;
                        local_msord++;
                        vvec += avec;
                    } else {
                        uvec -= alpha * vvec;
                    }
                }
                // If not solved, restart (outer loop continues)
            }

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
