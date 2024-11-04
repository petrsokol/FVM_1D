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

double computeTimeStep(double CFL, const std::vector<Conservative> &w, double *t, double dx);
double computeRezi(const std::vector<Conservative> &w, const std::vector<Conservative> &wn, double dt, double dx);
double bar(double rho_l, double rho_r, double vl, double vr);

Conservative flux(Conservative w, double q, double p);
Conservative fluxStar(Conservative w, double q, double S, double SM, double p, double p_star);

Conservative HLL(Conservative wl, Conservative wr);
Conservative HLLC(Conservative wl, Conservative wr);

Conservative minmod(Conservative a, Conservative b);
double minmod(double a, double b);

void exportData(const std::string &path, const std::string &filename, const std::vector<Conservative> &w, double dx);
std::string fileName(bool is_hll, bool second_order, int snapshotCount, int innerCells);
Conservative TVD(int i, const std::vector<Conservative> &w, short sign, double dx);
Conservative ENO(int i, const std::vector<Conservative> &w);

double setInitialConditions(std::vector<Conservative> &w, std::vector<Conservative> &wn, Conservative W_L, Conservative W_R);

void updateBounds(const std::vector<Conservative> &w, std::vector<Conservative> &wn);

void computeScheme(const std::vector<Conservative> &w, std::vector<Conservative> &wn,
                   double dt, bool isHLL, bool isSecOrd, double dx);

void exportSnapshot(double t, int *snapshotCount, const std::vector<Conservative> &w,
                    bool isHLL, bool isSecOrd, double dx);

void runExperiment(int innerCells, bool isHLL, bool isSecOrd);

int getGhostLayers();

double getKappa();

int main() {
    int innerCells = 500;
    runExperiment(innerCells, true, false);
    runExperiment(innerCells, false, false);
    runExperiment(innerCells, true, true);
    runExperiment(innerCells, false, true);

    return EXIT_SUCCESS;
    /*
     * todo
     *   - automatickej output tak, abych mohl sestavit grafy v různých časových intervalech
     *   - minmod a tvd rekonstrukce druhýho řádu přesnosti
     *   - přepínání na libovolnou konfiguraci schématu a řádu přesnosti
     */
}

void runExperiment(const int innerCells, const bool isHLL, const bool isSecOrd) {

    constexpr double T = 0.2;
    std::vector<Conservative> w(innerCells + 2 * getGhostLayers());
    std::vector<Conservative> wn(innerCells + 2 * getGhostLayers());

    // set initial condition
    Conservative W_L = Conservative(0.445, 0.311, 8.9280);
    Conservative W_R = Conservative(0.5, 0, 1.4275);
    double dx = setInitialConditions(w, wn, W_L, W_R);

    double t = 0;
    double rezi = 1;
    // int snapshotCount = 0;
    while(t < T) {
        double dt = computeTimeStep(0.5, w, &t, dx);

        computeScheme(w, wn, dt, isHLL, isSecOrd, dx);

        updateBounds(w, wn);

        rezi = computeRezi(w, wn, dt, dx);

        // update totalCells
        w = wn;

        // exportSnapshot(t, &snapshotCount, w, dx);
    }

    exportData(path, fileName(isHLL, isSecOrd, 20, innerCells), w, dx);
}

void computeScheme(const std::vector<Conservative> &w, std::vector<Conservative> &wn, const double dt, const bool isHLL,
                   const bool isSecOrd, const double dx)
{
    const int GHOST_LAYERS = getGhostLayers();
    const int firstInner = GHOST_LAYERS;
    const int innerCells = (int) w.size() - 2 * GHOST_LAYERS;
    for (int i = firstInner; i < firstInner + innerCells; ++i) {
        Conservative wl;
        Conservative wr;
        int l = i - 1;
        int r = i;
        if (isSecOrd) {
            wl = TVD(l, w, 1, dx);
            wr = TVD(r, w, -1, dx);
        } else {
            wl = w.at(l);
            wr = w.at(r);
        }
        Conservative flux = isHLL ? HLL(wl, wr) : HLLC(wl, wr);
        wn.at(l) -= dt / dx * flux;
        wn.at(r) += dt / dx * flux;
    }
}

int getGhostLayers() {
    return 2;
}

void updateBounds(const std::vector<Conservative> &w, std::vector<Conservative> &wn) {
    const int firstInner = getGhostLayers();
    const int innerCells = (int) w.size() - 2 * getGhostLayers();

    wn.at(0) = w.at(firstInner);
    wn.at(1) = w.at(firstInner);

    wn.at(firstInner + innerCells - 1) = w.at( firstInner + innerCells - 2); // first ghost cell = last inner cell
    wn.at(firstInner + innerCells - 2) = w.at(firstInner + innerCells - 2); // second ghost cell = last inner cell

}

double setInitialConditions(std::vector<Conservative> &w, std::vector<Conservative> &wn, const Conservative W_L, const Conservative W_R) {
    // returns value of dx;
    constexpr double X_LOWER_BOUND = 0;
    constexpr double X_UPPER_BOUND = 2;
    constexpr double X_MIDDLE = 0.7;

    constexpr int GHOST_LAYERS = 2;
    const int innerCells = (int) w.size() - 2 * GHOST_LAYERS;
    double dx = (X_UPPER_BOUND - X_LOWER_BOUND) / innerCells;

    int cells = innerCells + 2 * GHOST_LAYERS;

    for (int i = 0; i < cells; ++i) {
        double x = (i - GHOST_LAYERS) * dx;
        if (x < X_MIDDLE) {
            w.at(i) = W_L;
        } else {
            w.at(i) = W_R;
        }
    }

    wn = w;

    return dx;
}

Conservative ENO(int i, const std::vector<Conservative> &w) {
    Conservative res;
    Conservative a;
    return res;
}

Conservative TVD(int i, const std::vector<Conservative> &w, short sign, const double dx) {
    Conservative res;
    Conservative a = (w.at(i + 1) - w.at(i)) / dx;
    Conservative b = (w.at(i) - w.at(i - 1)) / dx;
    Conservative sigma_i = minmod(a, b);
    res = w.at(i) + sign * dx / 2 * sigma_i;
    return res;
}

double minmod(double a, double b) {
    if (a * b <= 0) {
        return 0;
    } else if (fabs(a) <= fabs(b)) {
        return a;
    } else if (fabs(a) > fabs(b)) {
        return b;
    } else {
        printf("sus: a: %f, b: %f", a, b);
        return -1;
    }
}

Conservative minmod(Conservative a, Conservative b) {
    Conservative res;
    res.r1 = minmod(a.r1, b.r1);
    res.r2 = minmod(a.r2, b.r2);
    res.r3 = minmod(a.r3, b.r3);
    return res;
}

std::string fileName(bool is_hll, bool second_order, int snapshotCount, int innerCells) {
    std::string scheme = is_hll ? "HLL" : "HLLC";
    std::string order = second_order == 1 ? "higherOrder" : "firstOrder";
    std::string slice = std::to_string(snapshotCount);
    std::string innerCellCount = std::to_string(innerCells);
    std::string name = order + "_" + scheme + "_" + slice + "_" + innerCellCount;
    return name;
}

double computeTimeStep(double CFL, const std::vector<Conservative> &w, double *t, const double dx) {
    double res = 1e9;
    Primitive pv{};
    const int cells = (int) w.size();
    for (int i = 0; i < cells; ++i) {
        pv = Primitive::computePV(w.at(i));
        res = fmin(res, dx / (fabs(pv.u * pv.c)));
    }
    *t += res;
    return res;
}

double computeRezi(const std::vector<Conservative> &w, const std::vector<Conservative> &wn, double dt, const double dx) {
    double res = 0;
    const int cells = (int) w.size();
    for (int i = 0; i < cells; ++i) {
        res += pow((wn.at(i).r1 - w.at(i).r1) / dt, 2) * dx;
    }
    return log10(sqrt(res));
}

double bar(double rho_l, double rho_r, double vl, double vr) {
    return (sqrt(rho_l) * vl + sqrt(rho_r) * vr) / (sqrt(rho_l) + sqrt(rho_r));
}

Conservative HLLC(Conservative wl, Conservative wr) {
    Conservative res;
    Primitive pvl = Primitive::computePV(wl);
    Primitive pvr = Primitive::computePV(wr);
    const double KAPPA = getKappa();

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
        res = flux( wl, ql, pvl.p);
    } else if (SL <= 0 && 0 < SM) {
        res = flux( wlStar, SM, p_star);
    } else if (SM <= 0 && 0 <= SR) {
        res = flux( wrStar, SM, p_star);
    } else if (SR < 0) {
        res = flux( wr, qr, pvr.p);
    } else {
        std::cout << "Scheme::HLLC: unreachable state \n";
    }

    return res;
}

double getKappa() {
    return 1.4;
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
    }

    return res;
}

Conservative fluxStar(Conservative w, double q, double S, double SM, double p, double p_star) {
    Conservative res{};

    res.r1 = w.r1 * (S - q) + 0;
    res.r2 = w.r2 * (S - q) + (p_star - p);
    res.r3 = w.r3 * (S - q) + p_star * SM - p * q;

    return res;
}

Conservative flux(Conservative w, double q, double p) {
    Conservative res{};

    res.r1 = w.r1 * q;
    res.r2 = w.r2 * q + p;
    res.r3 = (w.r3 + p) * q;

    return res;
}

void exportData(const std::string &path, const std::string &filename, const std::vector<Conservative> &w, const double dx) {
    std::ofstream output(path + "\\" + filename + ".dat");
    const int innerCells = (int) w.size() - 2 * getGhostLayers();
    for (int i = 0; i < innerCells; ++i) {
        output << i * dx << " " << w[i].r1 << " " << w[i].r2 << " " << w[i].r3 << std::endl;
    }
    output.close();
    std::cout << "Exported successfully" << std::endl;
}