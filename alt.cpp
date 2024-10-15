//
// Created by petrs on 15.10.2024.
//

#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include "structures/Conservative.h"
#include "structures/Primitive.h"

const std::string path = R"(C:\Users\petrs\Documents\CTU\BP\Charts\Data)";
constexpr int GHOST_LAYERS = 1;
constexpr int KAPPA = 1.4;

constexpr int innerCells = 150;
constexpr int cells = GHOST_LAYERS + innerCells + GHOST_LAYERS;
constexpr int firstInner = GHOST_LAYERS;
constexpr double T = 0.2;
constexpr double xLowerBound = 0;
constexpr double xUpperBound = 1.5;
constexpr double initialConditionBound = 0.7;
constexpr double dx = (xUpperBound - xLowerBound) / innerCells;


double computeTimeStep(double CFL, std::vector<Conservative> w);
double computeRezi(std::vector<Conservative> w, std::vector<Conservative> wn, double dt);
double bar(double rho_l, double rho_r, double vl, double vr);

Conservative flux(Conservative w, double q, double p);
Conservative fluxStar(Conservative w, double q, double S, double SM, double p, double p_star);
Conservative HLL(Conservative wl, Conservative wr);

Conservative HLLC(Conservative wl, Conservative wr);

void exportData(const std::string& path, const std::string& filename, std::vector<Conservative> w) {
    std::ofstream output(path + "\\" + filename + ".dat");
    for (int i = 0; i < innerCells; ++i) {
        output << w[i].r1 << std::endl;
    }
    output.close();
    std::cout << "Exported successfully" << std::endl;
}

int main() {
    std::vector<Conservative> w(cells);
    std::vector<Conservative> wn(cells);

    // set initial condition
    Conservative wl_initial = Conservative(0.445, 0.311, 8.9280);
    Conservative wr_initial = Conservative(0.5, 0, 1.4275);
    for (int i = 0; i < cells; ++i) {
        double x = (i - GHOST_LAYERS) * dx;
        if (x < initialConditionBound) {
            w.at(i) = wl_initial;
        } else {
            w.at(i) = wr_initial;
        }
    }
    wn = w;

    double t = 0;
    double rezi = 1;
    while(t < T) {

        // compute time step
        double dt = computeTimeStep(0.8, w);
        t += dt;

        // compute scheme
        for (int i = 0; i < cells - 1; ++i) {
            int l = i;
            int r = i + 1;
            Conservative wl = w.at(l);
            Conservative wr = w.at(r);
            Conservative flux = HLLC(wl, wr);

            wn.at(l) -= dt / dx * flux;
            wn.at(r) += dt / dx * flux;
        }

        // update bounds
        wn.at(0) = w.at(firstInner);
        wn.at(cells - 1) = w.at(cells - 2);

        // compute rezi
        rezi = computeRezi(w, wn, dt);

        // update cells
        w = wn;
    }

    exportData(path, "hll_1D_alt", w);
    std::cout << "final time: " << T << std::endl;

    return EXIT_SUCCESS;
    /*
     * todo
     *   - automatickej output tak, abych mohl sestavit grafy v různých časových intervalech
     *   - minmod a tvd rekonstrukce druhýho řádu přesnosti
     *   - přepínání na libovolnou konfiguraci schématu a řádu přesnosti
     */
}

Conservative HLLC(Conservative wl, Conservative wr) {
    Primitive pvl = Primitive::computePV(wl);
    Primitive pvr = Primitive::computePV(wr);

    double ql = pvl.u;
    double qr = pvr.u;

    double h_bar = bar(pvl.rho, pvr.rho, pvl.h, pvr.h);
    double u_bar = bar(pvl.rho, pvr.rho, pvl.u, pvr.u);
    double U_bar_sq = pow(u_bar, 2); // certified J. Holman verze
    double c_bar = sqrt((KAPPA - 1) * (h_bar - 0.5 * U_bar_sq));
    double q_bar = u_bar;

    double lambda_1 = ql - pvl.c;
    double lambda_m = qr + pvr.c;
    double lambda_1Roe = q_bar - c_bar;
    double lambda_mRoe = q_bar + c_bar;

    double SL = fmin(lambda_1, lambda_1Roe);
    double SR = fmax(lambda_m, lambda_mRoe);
    double SM = (pvr.rho * qr * (SR - qr) - pvl.rho * ql * (SL - ql) + pvl.p - pvr.p) / (pvr.rho * (SR - qr) - pvl.rho * (SL - ql));

    double p_star = pvl.rho * (ql - SL) * (ql - SM) + pvl.p;

    Conservative wlStar = 1 / (SL - SM) * fluxStar(wl, ql, SL, SM, pvl.p, p_star);
    Conservative wrStar = 1 / (SR - SM) * fluxStar(wr, qr, SR, SM, pvr.p, p_star);

    if (SL > 0) {
        return flux( wl, ql, pvl.p);
    } else if (SL <= 0 && 0 < SM) {
        return flux( wlStar, SM, p_star);
    } else if (SM <= 0 && 0 <= SR) {
        return flux( wrStar, SM, p_star);
    } else if (SR < 0) {
        return flux( wr, qr, pvr.p);
    } else {
        std::cout << "Scheme::HLLC: unreachable state \n";
    }
}

Conservative fluxStar(Conservative w, double q, double S, double SM, double p, double p_star) {
    Conservative res{};

    res.r1 = w.r1 * (S - q) + 0;
    res.r2 = w.r2 * (S - q) + (p_star - p);
    res.r3 = w.r3 * (S - q) + p_star * SM - p * q;
    return res;
}

double bar(double rho_l, double rho_r, double vl, double vr) {
    return (sqrt(rho_l) * vl + sqrt(rho_r) * vr) / (sqrt(rho_l) + sqrt(rho_r));
}

double computeTimeStep(double CFL, std::vector<Conservative> w) {
    double res = 1e9;
    Primitive pv{};
    for (int i = 0; i < cells; ++i) {
        pv = Primitive::computePV(w.at(i));
        res = fmin(res, dx / (fabs(pv.u * pv.c)));
    }
    return res;
}

double computeRezi(std::vector<Conservative> w, std::vector<Conservative> wn, double dt) {
    double res = 0;
    for (int i = 0; i < cells; ++i) {
        res += pow((wn.at(i).r1 - w.at(i).r1) / dt, 2) * dx;
    }
    return log10(sqrt(res));
}

Conservative HLL(Conservative wl, Conservative wr) {
    Conservative res{};

    Primitive pvl = Primitive::computePV(wl);
    Primitive pvr = Primitive::computePV(wr);

    double ql = pvl.u; // normálová rychlost
    double qr = pvr.u; // DP - \tilde u

    double SL = fmin(ql - pvl.c, qr - pvr.c);
    double SR = fmax(ql + pvl.c, qr + pvr.c);

    Conservative FL = flux(wl, ql, pvl.p);
    Conservative FR = flux(wr, qr, pvr.p);

    if (SL > 0) {
        res = FL;
    } else if (SL <= 0 && 0 <= SR) {
        res = (SR * FL - SL * FR + SR * SL * (wr - wl)) / (SR - SL);
    } else if (SR < 0) {
        res = FR;
    } else {
        printf("Error in HLL scheme, SL = %f, SR = %f\n", SL, SR);
        wl.toString();
        wr.toString();

    }

    return res;
}

Conservative flux(Conservative w, double q, double p) {
    Conservative res{};

    res.r1 = w.r1 * q;
    res.r2 = w.r2 * q + p;
    res.r3 = (w.r3 + p) * q;

    return res;
}