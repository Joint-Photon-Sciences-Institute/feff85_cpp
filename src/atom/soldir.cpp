// Dirac equation solver with iterative energy adjustment.
// Converted from: src/ATOM/soldir.f
//
// Solves the radial Dirac equation:
//   p' - kap*p/r = -(en/cl - v)*g - eg/r
//   g' + kap*g/r = (2*cl + en/cl - v)*p + ep/r
//
// INDEX CONVENTION: All loop indices (mat, max0, np, i, j, etc.) remain
// 1-based as in the Fortran original. Array accesses use the F() macro
// to convert to 0-based C++ indexing.

#include "soldir.hpp"
#include "intdir.hpp"
#include <cmath>
#include <cstring>

namespace feff::atom {

// 1-based to 0-based indexing
#define F(arr, i) (arr)[(i) - 1]

/// Calculate norm b. Extracted as separate subroutine in Fortran.
static void norm_calc(double& b, double hp[], const double dr[],
                      const double gg[], const double gp[],
                      const double ag[], const double ap[],
                      int method, double hx, int ndor,
                      double gpmat, double fl, int max0, int mat)
{
    b = 0.0;
    for (int i = 1; i <= max0; ++i) {
        F(hp, i) = F(dr, i) * (F(gg, i) * F(gg, i) + F(gp, i) * F(gp, i));
    }
    if (method == 1) {
        F(hp, mat) += F(dr, mat) * (gpmat * gpmat - F(gp, mat) * F(gp, mat)) / 2.0;
    }
    for (int i = 2; i <= max0; i += 2) {
        b += F(hp, i) + F(hp, i) + F(hp, i + 1);
    }
    b = hx * (b + b + F(hp, 1) - F(hp, max0)) / 3.0;
    for (int i = 1; i <= ndor; ++i) {
        double g = fl + fl + i;
        g = std::pow(F(dr, 1), g) / g;
        for (int j = 1; j <= i; ++j) {
            b += F(ag, j) * g * F(ag, i + 1 - j) + F(ap, j) * g * F(ap, i + 1 - j);
        }
    }
}

void soldir(double& en, double fl, double& agi, double& api, double& ainf,
            int nq, int kap, int& max0, int& ifail,
            DiracWorkspaceReal& work, MeshParamsReal& mesh,
            DiracSolverState& solver, ErrorState& error)
{
    static const char dprlab[9] = "  soldir";
    static const char drplab[9] = "  intdir";

    // Local working arrays
    double hg[atom_grid], agh[max_dev], hp_arr[atom_grid], aph[max_dev];
    double bg_arr[atom_grid], bgh[max_dev], bp_arr[atom_grid], bph[max_dev];

    // Aliases into workspace
    double* gg  = work.dg;
    double* ag  = work.ag;
    double* gp  = work.dp;
    double* ap  = work.ap;
    double* dv  = work.dv;
    double* av  = work.av;
    double* eg  = work.eg;
    double* ceg = work.ceg;
    double* ep  = work.ep;
    double* cep = work.cep;
    double* dr  = mesh.dr;
    int np   = mesh.np;
    int ndor = mesh.ndor;
    int nes  = mesh.nes;
    int idim = mesh.idim;

    double& ell  = solver.ell;
    double& fk   = solver.fk;
    double& ccl  = solver.ccl;
    int& imm     = solver.imm;
    int& nd      = solver.nd;
    int& node    = solver.node;
    int& mat     = solver.mat;

    std::memcpy(error.dlabpr, dprlab, 9);

    double enav = 1.0;
    ainf = std::abs(ainf);
    ccl = work.cl + work.cl;
    int iex = mesh.method;
    if (mesh.method <= 0) mesh.method = 1;

    fk = kap;
    if (av[0] < 0.0 && kap > 0) api = -agi * (fk + fl) / av[0];
    if (av[0] < 0.0 && kap < 0) api = -agi * av[0] / (fk - fl);
    ell = fk * (fk + 1.0) / ccl;
    node = nq - std::abs(kap);
    if (kap < 0) node = node + 1;

    double emin = 0.0;
    for (int i = 1; i <= np; ++i) {
        double a = (ell / (F(dr, i) * F(dr, i)) + F(dv, i)) * work.cl;
        if (a < emin) emin = a;
    }
    if (emin >= 0.0) {
        error.numerr = 75011;
        return;
    }
    if (en < emin) en = emin * 0.9;
    double edep = en;

    // Variables that persist across goto targets -- declared before any labels
    double b = 0.0, c = 0.0, f = 0.0, g = 0.0;
    double ggmat = 0.0, gpmat = 0.0, hgmat = 0.0, hpmat = 0.0;
    double bgmat = 0.0, bpmat = 0.0;
    double ah = 0.0;
    double test = 0.0;
    double einf = 0.0, esup = 0.0;
    int ies = 0;
    int j = 0;
    int jes = 0;
    int modmat = 0;
    double a_loc = 0.0;  // local 'a' used in long sections

label_101:
    error.numerr = 0;
    test = mesh.test1;
    if (mesh.method > 1) test = mesh.test2;
    einf = 1.0;
    esup = emin;
    en = edep;
    ies = 0;
    nd = 0;

label_105:
    jes = 0;

label_106:
    modmat = 0;
    imm = 0;
    if (std::abs((enav - en) / en) < 1.0e-01) imm = 1;
    enav = en;

label_107:
    // Integration of the inhomogeneous system
    for (int i = 1; i <= idim; ++i) {
        F(gg, i) = F(eg, i);
        F(gp, i) = F(ep, i);
    }
    for (int i = 2; i <= ndor; ++i) {
        F(ag, i) = F(ceg, i - 1);
        F(ap, i) = F(cep, i - 1);
    }
    intdir(gg, gp, ag, ap, ggmat, gpmat, en, fl, agi, api, ainf, max0,
           work, mesh, solver, error);
    if (error.numerr != 0) {
        std::memcpy(error.dlabpr, drplab, 9);
        return;
    }
    if (iex == 0) {
        // Match large component for homogeneous system (method=0)
        a_loc = ggmat / F(gg, mat);
        for (int i = mat; i <= max0; ++i) {
            F(gg, i) = a_loc * F(gg, i);
            F(gp, i) = a_loc * F(gp, i);
        }
        j = mat;
        goto label_215;
    }

    // Integration of the homogeneous system (label_141)
    for (int i = 1; i <= idim; ++i) {
        F(hg, i) = 0.0;
        F(hp_arr, i) = 0.0;
    }
    for (int i = 1; i <= ndor; ++i) {
        F(agh, i) = 0.0;
        F(aph, i) = 0.0;
    }
    imm = 1;
    if (mesh.method == 1) imm = -1;
    intdir(hg, hp_arr, agh, aph, hgmat, hpmat, en, fl, agi, api, ainf, max0,
           work, mesh, solver, error);

    // Match the large component
    a_loc = F(gg, mat) - ggmat;
    if (mesh.method < 2) {
        b = -a_loc / F(hg, mat);
    } else {
        b = F(gp, mat) - gpmat;
        ah = hpmat * F(hg, mat) - hgmat * F(hp_arr, mat);
        if (ah == 0.0) goto label_263;
        c = (b * F(hg, mat) - a_loc * F(hp_arr, mat)) / ah;
        b = (b * hgmat - a_loc * hpmat) / ah;
        for (int i = 1; i <= ndor; ++i) {
            F(ag, i) += c * F(agh, i);
            F(ap, i) += c * F(aph, i);
        }
        j = mat - 1;
        for (int i = 1; i <= j; ++i) {
            F(gg, i) += c * F(hg, i);
            F(gp, i) += c * F(hp_arr, i);
        }
    }
    for (int i = mat; i <= max0; ++i) {
        F(gg, i) += b * F(hg, i);
        F(gp, i) += b * F(hp_arr, i);
    }

    if (mesh.method >= 2) {
        // Integration of the system derived from disagreement in energy
        for (int i = 2; i <= ndor; ++i) {
            F(bgh, i) = F(ag, i - 1) / work.cl;
            F(bph, i) = F(ap, i - 1) / work.cl;
        }
        for (int i = 1; i <= max0; ++i) {
            F(bg_arr, i) = F(gg, i) * F(dr, i) / work.cl;
            F(bp_arr, i) = F(gp, i) * F(dr, i) / work.cl;
        }
        intdir(bg_arr, bp_arr, bgh, bph, bgmat, bpmat, en, fl, agi, api, ainf, max0,
               work, mesh, solver, error);

        // Match both components (method=2)
        f = F(bg_arr, mat) - bgmat;
        g = F(bp_arr, mat) - bpmat;
        a_loc = (g * F(hg, mat) - f * F(hp_arr, mat)) / ah;
        g = (g * hgmat - f * hpmat) / ah;
        for (int i = 1; i <= j; ++i) {
            F(bg_arr, i) += a_loc * F(hg, i);
            F(bp_arr, i) += a_loc * F(hp_arr, i);
        }
        for (int i = 1; i <= ndor; ++i) {
            F(bgh, i) += a_loc * F(agh, i);
            F(bph, i) += a_loc * F(aph, i);
        }
        for (int i = mat; i <= max0; ++i) {
            F(bg_arr, i) += g * F(hg, i);
            F(bp_arr, i) += g * F(hp_arr, i);
        }

        // Calculate the norm
        norm_calc(b, hp_arr, dr, gg, gp, ag, ap,
                  mesh.method, mesh.hx, ndor, gpmat, fl, max0, mat);

        // Correction to the energy (method=2)
        for (int i = 1; i <= max0; ++i) {
            F(hg, i) = (F(gg, i) * F(bg_arr, i) + F(gp, i) * F(bp_arr, i)) * F(dr, i);
        }
        ah = 0.0;
        c = 0.0;
        for (int i = 2; i <= max0; i += 2) {
            ah += F(hg, i) + F(hg, i) + F(hg, i + 1);
        }
        ah = mesh.hx * (ah + ah + F(hg, 1) - F(hg, max0)) / 3.0
             + F(hg, 1) / (fl + fl + 1.0);
        f = (1.0 - b) / (ah + ah);
        c = 1.0 - b;
        for (int i = 1; i <= max0; ++i) {
            F(gg, i) += f * F(bg_arr, i);
            F(gp, i) += f * F(bp_arr, i);
        }
        for (int i = 1; i <= ndor; ++i) {
            F(ag, i) += f * F(bgh, i);
            F(ap, i) += f * F(bph, i);
        }
    }

    // Search for the maximum of the modulus of large component
    a_loc = 0.0;
    F(bgh, 1) = b;
    F(bph, 1) = ah;
    for (int i = 1; i <= max0; ++i) {
        g = F(gg, i) * F(gg, i);
        if (g > a_loc) {
            a_loc = g;
            j = i;
        }
    }
    if (j > mat && modmat == 0) {
        modmat = 1;
        mat = j;
        if (mat % 2 == 0) mat = mat + 1;
        imm = 1;
        if (mat < (max0 - 10)) goto label_107;
        mat = max0 - 12;
        j = mat;
        if (mat % 2 == 0) mat = mat + 1;
    }

    // Compute number of nodes
label_215:
    nd = 1;
    j = (j > mat) ? j : mat;
    for (int i = 2; i <= j; ++i) {
        if (F(gg, i - 1) != 0.0) {
            if ((F(gg, i) / F(gg, i - 1)) <= 0.0) nd = nd + 1;
        }
    }

    if (nd == node) goto label_300;
    if (nd < node) {
        esup = en;
        if (einf < 0.0) goto label_271;
        en = en * 8.0e-01;
        if (std::abs(en) > mesh.test1) goto label_285;
        error.numerr = 238031;
        goto label_899;
    }
    einf = en;
    if (esup > emin) goto label_271;

label_263:
    en = en * 1.2;
    if (en > emin) goto label_285;
    error.numerr = 245041;
    goto label_899;

label_271:
    if (std::abs(einf - esup) <= mesh.test1) {
        error.numerr = 249051;
        goto label_899;
    }
    en = (einf + esup) / 2.0;

label_285:
    jes = jes + 1;
    if (jes <= nes) goto label_106;

    // Warning: too many attempts to find correct node count
    ifail = 1;

label_300:
    // Calculation of the norm
    norm_calc(b, hp_arr, dr, gg, gp, ag, ap,
              mesh.method, mesh.hx, ndor, gpmat, fl, max0, mat);

    if (mesh.method == 1) {
        // Correction to the energy (method=1)
        c = gpmat - F(gp, mat);
        f = F(gg, mat) * c * work.cl / b;
        if (gpmat != 0.0) c = c / gpmat;
    }

    en = en + f;
    g = std::abs(f / (en - f));

label_371:
    if ((en >= 0.0 || g > 2.0e-01) ||
        (std::abs(c) > test && (en < esup || en > einf))) {
        f = f / 2.0;
        g = g / 2.0;
        en = en - f;
        if (g > mesh.test1) goto label_371;
        error.numerr = 29071;
        goto label_899;
    }

    if (std::abs(c) > test) {
        ies = ies + 1;
        if (ies <= nes) goto label_105;
        ifail = 1;
    }

    // Divide by square root of the norm, and test sign of w.f.
    b = std::sqrt(b);
    c = b;
    if ((F(ag, 1) * agi) < 0.0 || (F(ap, 1) * api) < 0.0) c = -c;
    for (int i = 1; i <= ndor; ++i) {
        F(ag, i) /= c;
        F(ap, i) /= c;
    }
    if ((F(gg, 1) * agi) < 0.0 || (F(gp, 1) * api) < 0.0) b = -b;
    for (int i = 1; i <= max0; ++i) {
        F(gg, i) /= b;
        F(gp, i) /= b;
    }
    if (max0 >= np) return;
    for (int i = max0 + 1; i <= np; ++i) {
        F(gg, i) = 0.0;
        F(gp, i) = 0.0;
    }
    return;

    // Abnormal exit
label_899:
    if (iex == 0 || mesh.method == 2) return;
    mesh.method = mesh.method + 1;
    goto label_101;
}

#undef F

} // namespace feff::atom
