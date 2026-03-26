// Dirac equation integrator using 5-point predictor-corrector method.
// Converted from: src/ATOM/intdir.f
//
// Numerical method:
//   predicted value   p(n) = y(n-1) + c * sum_{i=1,5} cop(i)*y'(n-i)
//   corrected value   c(n) = y(n-1) + c * sum_{i=1,4} coc(i)*y'(n-i) + coc(5)*p'(n)
//   final value       y(n) = cmix*c(n) + (1-cmix)*p(n)
//                     cmix = cmixn/cmixd
//
// INDEX CONVENTION: All internal loop indices (mat, max0, np, i, j, k)
// remain 1-based as in the Fortran original. Array accesses use the F()
// macro to convert to 0-based C++ indexing.

#include "intdir.hpp"
#include <cmath>

namespace feff::atom {

// 1-based to 0-based index helper
#define F(arr, idx) (arr)[(idx) - 1]

void intdir(double gg[], double gp[], double ag[], double ap[],
            double& ggmat, double& gpmat,
            double en, double fl, double agi, double api, double ainf,
            int& max0,
            DiracWorkspaceReal& work, MeshParamsReal& mesh,
            DiracSolverState& solver, ErrorState& error)
{
    // Predictor-corrector coefficients (Fortran DATA + SAVE)
    static double cop[5] = {2.51e+02, -1.274e+03, 2.616e+03, -2.774e+03, 1.901e+03};
    static double coc[5] = {-1.9e+01, 1.06e+02, -2.64e+02, 6.46e+02, 2.51e+02};
    static constexpr double cmixn = 4.73e+02;
    static constexpr double cmixd = 5.02e+02;
    static constexpr double hxd = 7.2e+02;
    static constexpr int npi = 5;
    static int icall = 0;
    static double cmc = 0.0;

    // Working arrays for the 5-point stencil derivatives
    double dg_w[5], dp_w[5];

    // Aliases
    double* dv = work.dv;
    double* av = work.av;
    double* dr = mesh.dr;
    int ndor = mesh.ndor;
    int np = mesh.np;
    double& fk  = solver.fk;
    double& ccl = solver.ccl;
    int& imm    = solver.imm;
    int& mat    = solver.mat;

    // First-call initialization: merge predictor and corrector coefficients
    if (icall == 0) {
        icall = 1;
        double cmix = cmixn / cmixd;
        double amix = 1.0 - cmix;
        // cmc = c * coc(5)  where c = cmixn/cmixd and coc(5) = 2.51e+02 initially
        cmc = cmix * coc[4];
        double f_val = coc[0];
        for (int j = 2; j <= npi; ++j) {
            double g_val = coc[j - 1];
            coc[j - 1] = cmix * f_val + amix * cop[j - 1];
            f_val = g_val;
        }
        coc[0] = cmix * cop[0];
    }

    double c = mesh.hx / hxd;
    double ec = en / work.cl;

    ag[0] = agi;
    ap[0] = api;

    // Variables for the integration loop -- declared here so goto doesn't
    // bypass initialization
    int i_pt = npi;  // current integration point (1-based)
    int k_dir = 1;   // direction: +1 outward, -1 inward
    double cmcc = 0.0;

    if (imm < 0) goto label_81;

    if (imm == 0) {
        // Search for the second sign-change point (classical turning point)
        mat = npi;
        int j_sign = 1;

    label_16:
        mat = mat + 2;
        if (mat >= np) {
            if (ec > -0.0003) {
                mat = np - 12;
                goto label_25;
            }
            error.numerr = 56011;
            return;
        }
        {
            double f_val = F(dv, mat) + solver.ell / (F(dr, mat) * F(dr, mat));
            f_val = (f_val - ec) * j_sign;
            if (f_val > 0.0) goto label_16;
        }

    label_25:
        j_sign = -j_sign;
        if (j_sign < 0) goto label_16;
        if (mat >= np - npi) mat = np - 12;
    }

    // Compute development coefficients at the origin
    for (int j = 2; j <= ndor; ++j) {
        int kk = j - 1;
        double a = fl + fk + kk;
        double b = fl - fk + kk;
        double ep_val = a * b + F(av, 1) * F(av, 1);
        double f_val = (ec + ccl) * F(ap, kk) + F(ap, j);
        double g_val = ec * F(ag, kk) + F(ag, j);
        for (int ii = 1; ii <= kk; ++ii) {
            f_val -= F(av, ii + 1) * F(ap, j - ii);
            g_val -= F(av, ii + 1) * F(ag, j - ii);
        }
        F(ag, j) = (b * f_val + F(av, 1) * g_val) / ep_val;
        F(ap, j) = (F(av, 1) * f_val - a * g_val) / ep_val;
    }

    // Initial values for outward integration at first npi points
    for (int ii = 1; ii <= npi; ++ii) {
        F(gg, ii) = 0.0;
        F(gp, ii) = 0.0;
        dg_w[ii - 1] = 0.0;
        dp_w[ii - 1] = 0.0;
        for (int j = 1; j <= ndor; ++j) {
            double a = fl + j - 1;
            double b = std::pow(F(dr, ii), a);
            double ab = a * b * c;
            F(gg, ii) += b * F(ag, j);
            F(gp, ii) += b * F(ap, j);
            dg_w[ii - 1] += ab * F(ag, j);
            dp_w[ii - 1] += ab * F(ap, j);
        }
    }

    i_pt = npi;
    k_dir = 1;
    ggmat = F(gg, mat);
    gpmat = F(gp, mat);

    // Integration of the inhomogeneous system
label_51:
    cmcc = cmc * c;

label_55:
    {
        double a = F(gg, i_pt) + dg_w[0] * cop[0];
        double b_val = F(gp, i_pt) + dp_w[0] * cop[0];
        i_pt = i_pt + k_dir;
        double ep_val = F(gp, i_pt);
        double eg_val = F(gg, i_pt);
        F(gg, i_pt) = a - dg_w[0] * coc[0];
        F(gp, i_pt) = b_val - dp_w[0] * coc[0];
        for (int j = 2; j <= npi; ++j) {
            a += dg_w[j - 1] * cop[j - 1];
            b_val += dp_w[j - 1] * cop[j - 1];
            F(gg, i_pt) += dg_w[j - 1] * coc[j - 1];
            F(gp, i_pt) += dp_w[j - 1] * coc[j - 1];
            dg_w[j - 2] = dg_w[j - 1];
            dp_w[j - 2] = dp_w[j - 1];
        }
        double f_val = (ec - F(dv, i_pt)) * F(dr, i_pt);
        double g_val = f_val + ccl * F(dr, i_pt);
        F(gg, i_pt) += cmcc * (g_val * b_val - fk * a + ep_val);
        F(gp, i_pt) += cmcc * (fk * b_val - f_val * a - eg_val);
        dg_w[npi - 1] = c * (g_val * F(gp, i_pt) - fk * F(gg, i_pt) + ep_val);
        dp_w[npi - 1] = c * (fk * F(gp, i_pt) - f_val * F(gg, i_pt) - eg_val);
    }
    if (i_pt != mat) goto label_55;

    if (k_dir < 0) return;  // label 999

    // Swap values at matching point
    {
        double a = ggmat;
        ggmat = F(gg, mat);
        F(gg, mat) = a;
        a = gpmat;
        gpmat = F(gp, mat);
        F(gp, mat) = a;
    }
    if (imm != 0) goto label_81;

    // Initial values for inward integration
    {
        double a = mesh.test1 * std::abs(ggmat);
        if (ainf > a) ainf = a;
    }
    max0 = np + 2;

label_73:
    {
        double a_limit = 7.0e+02 / work.cl;

    label_75:
        max0 = max0 - 2;
        if ((max0 + 1) <= (mat + npi)) {
            error.numerr = 138021;
            return;
        }
        if ((F(dv, max0) - ec) * F(dr, max0) * F(dr, max0) > a_limit) goto label_75;
    }

label_81:
    c = -c;
    {
        double a = -std::sqrt(-ec * (ccl + ec));
        if (a * F(dr, max0) < -1.7e+02) goto label_73;
        double b_val = a / (ccl + ec);
        double f_val = ainf / std::exp(a * F(dr, max0));
        if (f_val == 0.0) f_val = 1.0;
        for (int ii = 1; ii <= npi; ++ii) {
            int jj = max0 + 1 - ii;
            F(gg, jj) = f_val * std::exp(a * F(dr, jj));
            F(gp, jj) = b_val * F(gg, jj);
            dg_w[ii - 1] = a * F(dr, jj) * F(gg, jj) * c;
            dp_w[ii - 1] = b_val * dg_w[ii - 1];
        }
    }
    i_pt = max0 - npi + 1;
    k_dir = -1;
    goto label_51;
}

#undef F

} // namespace feff::atom
