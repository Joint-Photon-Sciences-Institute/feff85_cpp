// Graves-Morris Lanczos solver for FMS.
// Converted from: gggm.f
// NOTE: gggm uses A = T*G0 (positive sign), solving x - A*x = b via
// Graves-Morris Lanczos acceleration.
// Reference: Graves-Morris & Salam, Num.Algor. 21, p.213 (1999)

#include "gggm.hpp"

#include <cmath>
#include <complex>

namespace feff::fms {

static Complexf cdot_gm(int n, const Complexf* bra, const Complexf* ket) {
    Complexf sum(0.0f, 0.0f);
    for (int i = 0; i < n; ++i)
        sum += std::conj(bra[i]) * ket[i];
    return sum;
}

static bool converged_gm(int n, const Complexf* v, float tol) {
    for (int i = 0; i < n; ++i) {
        if (std::abs(v[i].real()) > tol || std::abs(v[i].imag()) > tol)
            return false;
    }
    return true;
}

void gggm(int nsp, const int* i0, int ipi, int ipf, const int* lipotx,
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

    // Build g0t = T*G0 (positive sign)
    g0t = Eigen::MatrixXcf::Zero(istate, istate);
    for (int icol = 0; icol < istate; ++icol) {
        for (int irow = 0; irow < istate; ++irow) {
            if (std::abs(g0(irow, icol)) > toler2) {
                g0t(irow, icol) += tmatrx_diag(irow) * g0(irow, icol);
            }
            int l1 = basis.l(irow);
            int m1 = basis.m(irow);
            int isp1 = basis.isp(irow);
            int m2 = m1 + isp1;
            if (nsp == 2 && m2 > -l1 + 1 && m2 < l1 + 2) {
                int ist2 = irow + ((isp1 == 0) ? 1 : -1);
                if (ist2 >= 0 && ist2 < istate &&
                    std::abs(g0(ist2, icol)) > toler2) {
                    g0t(irow, icol) += tmatrx_off(ist2) * g0(ist2, icol);
                }
            }
        }
    }

    Eigen::VectorXcf xvec(istate), wvec(istate), x0(istate), x1(istate);
    Eigen::VectorXcf bvec(istate), avec(istate);
    Eigen::VectorXcf r0(istate), r1(istate), t0(istate), t1(istate);

    for (int ip = ipi; ip <= ipf; ++ip) {
        int ipart = nsp * (lipotx[ip] + 1) * (lipotx[ip] + 1);
        for (int is1 = 0; is1 < ipart; ++is1) {
            int is2 = is1 + i0[ip];
            int l1_check = basis.l(is2);
            if (!lcalc[l1_check]) continue;

            xvec.setZero();
            bvec.setZero();
            bvec(is2) = Complexf(1.0f, 0.0f);
            int local_msord = 0;
            bool solved = false;
            Complexf q0_val(1.0f, 0.0f);

            // Restart loop
            for (int istart = 0; !solved; ++istart) {
                if (istart > 0) {
                    xvec += x0 / q0_val;
                    avec.noalias() = g0t.topLeftCorner(istate, istate) * xvec;
                    local_msord++;
                    bvec = avec - xvec;
                    bvec(is2) += Complexf(1.0f, 0.0f);
                }

                r0 = bvec;
                x0.setZero();
                x1 = bvec;

                r1.noalias() = g0t.topLeftCorner(istate, istate) * bvec;
                local_msord++;

                // Choose wvec for conditioning
                {
                    Complexf ww = cdot_gm(istate, r0.data(), r0.data());
                    Complexf aa = cdot_gm(istate, r1.data(), r1.data());
                    Complexf wa = cdot_gm(istate, r0.data(), r1.data());
                    Complexf aw = std::conj(wa);
                    Complexf dd = aa * ww - aw * wa;
                    if (std::abs(dd / (aa * ww)) < 1.0e-8f) {
                        wvec = r0 / ww;
                    } else {
                        Complexf ww2 = (ww - aw) / dd;
                        Complexf aa2 = (wa - aa) / dd;
                        wvec = r0 * aa2 + r1 * ww2;
                    }
                }

                Complexf e0 = cdot_gm(istate, wvec.data(), r0.data());
                Complexf e1 = cdot_gm(istate, wvec.data(), r1.data());
                q0_val = Complexf(1.0f, 0.0f);
                Complexf q1(1.0f, 0.0f);

                constexpr int nitx = 10;
                for (int nit = 1; nit <= nitx; ++nit) {
                    float tol = toler1 * std::abs(q1) / 10.0f;
                    if (converged_gm(istate, r1.data(), tol)) {
                        xvec += x1 / q1;
                        solved = true;
                        break;
                    }

                    Complexf alpha = e1 / e0;
                    t0 = r1 - alpha * r0;

                    t1.noalias() = g0t.topLeftCorner(istate, istate) * t0;
                    local_msord++;

                    Complexf wa2 = cdot_gm(istate, t0.data(), t1.data());
                    Complexf ww2 = cdot_gm(istate, t0.data(), t0.data());
                    Complexf aa2 = cdot_gm(istate, t1.data(), t1.data());
                    Complexf aw2 = std::conj(wa2);
                    Complexf theta = (wa2 - aa2) / (ww2 - aw2);

                    r0 = t1 - theta * t0;
                    Complexf dd = Complexf(1.0f, 0.0f) - theta;
                    x0 = t0 + dd * (x1 - alpha * x0);
                    q0_val = dd * (q1 - alpha * q0_val);

                    tol = toler1 * std::abs(q0_val);
                    if (converged_gm(istate, r0.data(), tol)) {
                        xvec += x0 / q0_val;
                        solved = true;
                        break;
                    }

                    // Prepare next iteration
                    e0 = cdot_gm(istate, wvec.data(), r0.data());
                    Complexf beta = e0 / e1;
                    t0 = r0 - beta * r1;

                    avec.noalias() = g0t.topLeftCorner(istate, istate) * t0;
                    local_msord++;

                    dd = beta * theta;
                    r1 = avec + dd * r1;
                    e1 = cdot_gm(istate, wvec.data(), r1.data());

                    dd = beta * (Complexf(1.0f, 0.0f) - theta);
                    x1 = x0 - dd * x1 + t0;
                    q1 = q0_val - (Complexf(1.0f, 0.0f) - theta) * beta * q1;
                }
                // If inner loop exhausted, restart (outer loop continues)
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
