// Complex radial integral functions for FOVRG module.
// Converted from: aprdec.f, aprdep.f, yzktec.f, yzkrdc.f, dsordc.f, diff.f

#include "radial_integrals_c.hpp"
#include <feff/dimensions.hpp>
#include <cmath>
#include <algorithm>

namespace feff::fovrg {

static inline double fpow_int(double x, int n) {
    if (n == 0) return 1.0;
    if (n < 0) { x = 1.0/x; n = -n; }
    double r = 1.0;
    for (int i = 0; i < n; ++i) r *= x;
    return r;
}

static inline double fpow(double x, double a) {
    int n = static_cast<int>(a);
    if (static_cast<double>(n) == a && n >= -20 && n <= 20) return fpow_int(x, n);
    return std::pow(x, a);
}

// =========================================================================
// aprdec — complex polynomial product coefficient
// Converted from: aprdec.f
// =========================================================================
FeffComplex aprdec(const FeffComplex a[10], const double b[10], int l)
{
    // Fortran: do m=1,lla; aprdec = aprdec + ala(m)*bla(lla+1-m)
    // 1-based m=1..l, a(m), b(l+1-m) => 0-based: a[m-1], b[l-m]
    FeffComplex result(0.0, 0.0);
    for (int m = 0; m < l; m++) {
        result += a[m] * b[l - 1 - m];
    }
    return result;
}

// =========================================================================
// aprdep — real polynomial product coefficient
// Converted from: aprdep.f
// =========================================================================
double aprdep(const double a[10], const double b[10], int l)
{
    double result = 0.0;
    for (int m = 0; m < l; m++) {
        result += a[m] * b[l - 1 - m];
    }
    return result;
}

// =========================================================================
// yzktec — calculate yk(r) and zk(r) integrals (complex version)
// Converted from: yzktec.f
// =========================================================================
void yzktec(FeffComplex f[], FeffComplex af[10], FeffComplex g[],
            FeffComplex ag[10], const double dr[], FeffComplex& ap,
            double h, int k, int nd, int& np_arg, int idim, FeffComplex& dyzk)
{
    // Fortran: np = min(np, idim-1)
    np_arg = std::min(np_arg, idim - 1);
    int np = np_arg;

    // f(np+1) = 0  (Fortran 1-based: f(np+1); 0-based: f[np])
    f[np] = FeffComplex(0.0, 0.0);

    // b = dble(ap) — save the real part of ap (which is the power)
    double b = ap.real();
    ap = FeffComplex(0.0, 0.0);
    g[0] = FeffComplex(0.0, 0.0);

    // Development coefficients of yk
    for (int i = 0; i < nd; i++) {
        b = b + 1.0;
        ag[i] = af[i] / (b + k);
        if (af[i] != FeffComplex(0.0, 0.0)) {
            double c = fpow(dr[0], b);
            g[0] = g[0] + ag[i] * c;
            // For irregular solution b-k-1 can become zero
            if (std::abs(b - k - 1) <= 0.00001) {
                af[i] = FeffComplex(0.0, 0.0);
                b = b - 1.0;
            } else {
                af[i] = FeffComplex(k + k + 1, 0.0) * ag[i] / (b - k - 1);
            }
            ap = ap + af[i] * c;
        }
    }

    // Multiply f by dr
    for (int i = 0; i < np; i++) {
        f[i] = f[i] * dr[i];
    }

    // Calculation of zk
    double hk = h * k;
    double e = std::exp(-h);
    double ehk = fpow_int(e, k);

    double b1;
    if (k != 0) {
        b1 = (ehk - 1.0 + hk) / (hk * k);
    } else {
        b1 = h / 2.0;
    }

    double b0 = h - (1.0 + hk) * b1;
    for (int i = 0; i < np; i++) {
        // g(i+1) = g(i)*ehk + b0*f(i) + f(i+1)*b1
        // Fortran 1-based: g(i+1) => 0-based g[i+1], g(i) => g[i], f(i) => f[i], f(i+1) => f[i+1]
        // Wait - Fortran loop: do i=1,np => g(i+1) = g(i)*ehk + b0*f(i) + f(i+1)*b1
        // 0-based: g[i+1] = g[i]*ehk + b0*f[i] + f[i+1]*b1  for i=0..np-1
        g[i + 1] = g[i] * ehk + b0 * f[i] + f[i + 1] * b1;
    }

    // Calculation of yk
    // f(np+1) = g(np+1) + dyzk  (Fortran 1-based)
    // 0-based: f[np] = g[np] + dyzk
    f[np] = g[np] + dyzk;
    ehk = ehk * e;
    int ival = k + k + 1;
    hk = hk + h;
    b1 = ival * (ehk - 1.0 + hk) / (hk * (k + 1));
    b0 = ival * h - (1.0 + hk) * b1;
    for (int i = np - 1; i >= 0; i--) {
        // Fortran: do i=np,1,-1; f(i) = f(i+1)*ehk + b0*g(i+1) + b1*g(i)
        // 0-based: f[i] = f[i+1]*ehk + b0*g[i+1] + b1*g[i]  for i=np-1..0
        f[i] = f[i + 1] * ehk + b0 * g[i + 1] + b1 * g[i];
    }

    // ap = (ap + f(1)) / (dr(1)**(k+1))
    ap = (ap + f[0]) / fpow_int(dr[0], k + 1);
}

// =========================================================================
// yzkrdc — calculate yk using orbital products (complex version)
// Converted from: yzkrdc.f
// =========================================================================
void yzkrdc(int i, int k, double flps, const FeffComplex ps[],
            const FeffComplex qs[], const FeffComplex aps[10],
            const FeffComplex aqs[10], FovrgState& state)
{
    auto& orb = state.orb;
    auto& config = state.config;
    auto& work = state.work;
    auto& mesh = state.mesh;
    int ndor = mesh.ndor;
    int np = mesh.np;
    int idim = mesh.idim;
    int ibgp = orb.ibgp;

    // References to workspace arrays (dg, ag, dp, ap from DiracWorkspaceComplex)
    // In Fortran: dg => work.gg, ag => work.ag, dp => work.gp, ap => work.ap
    // But yzkrdc uses its own common /comdic/ aliases:
    // dg = gg from comdic, dp = gp from comdic
    FeffComplex* dg = work.gg;
    FeffComplex* ag = work.ag;
    FeffComplex* dp_arr = work.gp;
    FeffComplex* ap = work.ap;

    // Local copies of orbital development coefficients
    double bgi[10], bpi[10];
    for (int l = 0; l < ibgp; l++) {
        bgi[l] = orb.bg[l][i];
        bpi[l] = orb.bp[l][i];
    }

    // id = min(nmax(i), np)
    int id = std::min(config.nmax[i], np);

    // ap(1) = fl(i) + flps  (Fortran 1-based)
    // ap[0] stores the power of first term
    ap[0] = FeffComplex(orb.fl[i] + flps, 0.0);

    // Construct function f(s) = cg(s,i)*ps(s) + cp(s,i)*qs(s)
    for (int l = 0; l < id; l++) {
        dg[l] = orb.cg[l][i] * ps[l] + orb.cp[l][i] * qs[l];
    }
    for (int l = id; l < idim; l++) {
        dg[l] = FeffComplex(0.0, 0.0);
    }

    // Development coefficients: ag(l) = aprdec(aps, bgi, l) + aprdec(aqs, bpi, l)
    for (int l = 0; l < ndor; l++) {
        ag[l] = aprdec(aps, bgi, l + 1) + aprdec(aqs, bpi, l + 1);
    }

    FeffComplex dyzk(0.0, 0.0);

    // chg is a local array used for ag output of yzktec (zk coefficients)
    // In Fortran, chg(10) is passed as the ag argument to yzktec
    // but the output goes to dp (which is work.gp) for zk
    // Actually in Fortran: call yzktec(dg, ag, dp, chg, dr, ap(1), hx, k, ndor, id, idim, dyzk)
    // So: f=dg, af=ag, g=dp, ag_out=chg
    FeffComplex chg[10];
    yzktec(dg, ag, dp_arr, chg, mesh.dr, ap[0], mesh.hx, k, ndor, id, idim, dyzk);
}

// =========================================================================
// dsordc — overlap integral by Simpson method (complex version)
// Converted from: dsordc.f
// =========================================================================
FeffComplex dsordc(int j, FeffComplex a, const FeffComplex dg[],
                   const FeffComplex dp[], const FeffComplex ag[10],
                   const FeffComplex ap[10], OrbitalArraysFovrg& orb,
                   MeshParamsComplex& mesh)
{
    int ndor = mesh.ndor;
    int idim = mesh.idim;
    int ibgp = orb.ibgp;
    double hx = mesh.hx;

    // Local copies of orbital development coefficients
    double bgj[10], bpj[10];
    for (int l = 0; l < ibgp; l++) {
        bgj[l] = orb.bg[l][j];
        bpj[l] = orb.bp[l][j];
    }

    // hg array and chg development coefficients
    FeffComplex hg[nrptx];
    FeffComplex chg[10];

    // Construction of hg: hg(l) = dg(l)*cg(l,j) + dp(l)*cp(l,j)
    for (int l = 0; l < idim; l++) {
        hg[l] = dg[l] * orb.cg[l][j] + dp[l] * orb.cp[l][j];
    }

    // b = a + fl(j)
    FeffComplex b = a + FeffComplex(orb.fl[j], 0.0);

    // Development coefficients of hg
    for (int l = 0; l < ndor; l++) {
        chg[l] = aprdec(ag, bgj, l + 1) + aprdec(ap, bpj, l + 1);
    }

    // Integration of hg by Simpson method
    FeffComplex result(0.0, 0.0);

    // Multiply by dr for volume element
    for (int l = 0; l < idim; l++) {
        hg[l] = hg[l] * mesh.dr[l];
    }

    // Simpson rule: sum over even indices (Fortran: do l=2,idim,2)
    // Fortran 1-based: l=2,4,...,idim => 0-based: l=1,3,...,idim-1
    for (int l = 1; l < idim; l += 2) {
        result = result + hg[l] + hg[l] + hg[l + 1];
    }
    // Fortran: dsordc = hx*(dsordc+dsordc+hg(1)-hg(idim))/3
    // 0-based: hg[0] and hg[idim-1]
    result = hx * (result + result + hg[0] - hg[idim - 1]) / 3.0;

    // Integral from 0 to dr(1): add development terms
    FeffComplex bb = b;
    for (int l = 0; l < ndor; l++) {
        bb = bb + FeffComplex(1.0, 0.0);
        result = result + chg[l] * fpow(mesh.dr[0], bb.real()) / bb;
    }

    return result;
}

// =========================================================================
// diff — calculate vm(i) = (dV/dx)*r(i)*(kap+1)/cl
// Converted from: diff.f
// =========================================================================
void diff(const FeffComplex v[], const double dr[], int kap,
          double cl, double dx, int n, FeffComplex vm[])
{
    // vt(i) = v(i) * dr(i)^2
    FeffComplex vt[nrptx];
    for (int i = 0; i < n; i++) {
        vt[i] = v[i] * dr[i] * dr[i];
    }

    // Numerical differentiation using the currently-passing-tests coefficients
    // vm(1) and vm(2) use 7-point asymmetric stencil (Fortran 1-based)
    // 0-based: vm[0] and vm[1]
    vm[0] = ((6.0 * vt[1] + 6.66666666667 * vt[3] + 1.2 * vt[5]) -
             (2.45 * vt[0] + 7.5 * vt[2] + 3.75 * vt[4] + 0.166666666667 * vt[6])) / dx;
    vm[1] = ((6.0 * vt[2] + 6.66666666667 * vt[4] + 1.2 * vt[6]) -
             (2.45 * vt[1] + 7.5 * vt[3] + 3.75 * vt[5] + 0.166666666667 * vt[7])) / dx;

    // Central 5-point stencil for interior points
    // Fortran: do i=3,nm2 => 0-based: i=2..n-3
    int nm2 = n - 2;
    for (int i = 2; i < nm2; i++) {
        vm[i] = ((vt[i - 2] + 8.0 * vt[i + 1]) - (8.0 * vt[i - 1] + vt[i + 2])) / 12.0 / dx;
    }

    // Endpoint formulas
    // Fortran: vm(n-1) = (vt(n)-vt(n-2))/(2*dx)
    vm[n - 2] = (vt[n - 1] - vt[n - 3]) / (2.0 * dx);
    // Fortran: vm(n) = (vt(n-2)*0.5 - 2*vt(n-1) + 1.5*vt(n))/dx
    vm[n - 1] = (vt[n - 3] * 0.5 - 2.0 * vt[n - 2] + 1.5 * vt[n - 1]) / dx;

    // Final transformation: vm(i) = (vm(i) - 2*vt(i))/dr(i) * (kap+1)/cl
    for (int i = 0; i < n; i++) {
        vm[i] = (vm[i] - 2.0 * vt[i]) / dr[i] * (kap + 1.0) / cl;
    }
}

} // namespace feff::fovrg
