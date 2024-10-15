//
// Created by petrs on 02.12.2023.
//

#include <cmath>
#include "Scheme.h"
#include "../structures/Primitive.h"
#include "Def.h"
#include "Bound.h"
#include "NACA.h"

void Scheme::updateCellDT(std::unordered_map<int, Cell> &cells, double CFL, bool useGlobalTimeStep) {

    double globalDT = 1e9;

    for (auto &cell: cells) {
        Primitive pv = Primitive::computePV(cell.second.w);

        double u_xi = fabs(pv.u * cell.second.xi.ux + pv.v * cell.second.xi.uy);
        double u_eta = fabs(pv.u * cell.second.eta.ux + pv.v * cell.second.eta.uy);

        double d_xi = (u_xi + pv.c) / cell.second.xi.length;
        double d_eta = (u_eta + pv.c) / cell.second.eta.length;

        double res = CFL / (d_xi + d_eta);

        if(useGlobalTimeStep) {
            globalDT = fmin(globalDT, res);
        } else {
            // local time step is assigned to each cell
            cell.second.dt = res;
        }
    }

    if(useGlobalTimeStep) {
        // global time step is assigned to every cell
        for (auto &cell: cells) {
            cell.second.dt = globalDT;
        }
    }
}

Conservative Scheme::HLL(const std::unordered_map<int, Cell> &cells, const Interface &face) {
    Conservative res{};
    Conservative wl = cells.at(face.left).w;
    Conservative wr = cells.at(face.right).w;
    Primitive pvl = Primitive::computePV(wl);
    Primitive pvr = Primitive::computePV(wr);

    if(_isnan(pvr.p)) {
        std::cout << "tlak je divnej";
        wr.toString();
        face.toString();
    }

    double ql = pvl.u * face.nx + pvl.v * face.ny; // normálová rychlost
    double qr = pvr.u * face.nx + pvr.v * face.ny;

    double SL = fmin(ql - pvl.c, qr - pvr.c);
    double SR = fmax(ql + pvl.c, qr + pvr.c);

    Conservative FL = Scheme::flux(face, wl, ql, pvl.p);
    Conservative FR = Scheme::flux(face, wr, qr, pvr.p);

    if (_isnan(FL.r1)) {
        std::cout << "error at left side " << face.left << std::endl;
        wl.toString();
    }
    if (_isnan(FR.r1)) {
//        std::cout << "error at right side " << face.right << ", is inner: " << cells.at(face.right).isInner << std::endl;
        Def::coordsToString(face.right);
//        wr.toString();
        Def::errorCount++;
    }

    if (SL > 0) {
        res = FL;
    } else if (SL <= 0 && 0 <= SR) {
        res = (SR * FL - SL * FR + SR * SL * (wr - wl)) / (SR - SL);
    } else if (SR < 0) {
        res = FR;
    } else {
        std::cout << "error in HLL scheme" << std::endl;
        Def::error = true;
    }

    return res;
}

void Scheme::computeScheme(std::unordered_map<int, Cell> &cells,
                           const std::unordered_map<std::pair<int, int>, Interface, pair_hash> &faces) {
    // loops over all interfaces
    for (const auto &face: faces) {
        Conservative flux = Def::isHLLC ? HLLC(cells, face.second) : HLL(cells, face.second);

        // selects cells left and right from the current interfaces, adds or subtracts flux
        cells.at(face.second.left).rezi -= cells.at(face.second.left).dt / cells.at(face.second.left).area * flux * face.second.length;
        cells.at(face.second.right).rezi += cells.at(face.second.right).dt / cells.at(face.second.right).area * flux * face.second.length;
    }

    bool test = true;
    if (test && Def::isNaca) {
        for (int i = 0; i < NACA::wingLength; ++i) {
            int k = Def::firstInner + NACA::wingStart + i;
            // face is defined by two end points, not the neighbouring cells
            Interface face = faces.at(std::make_pair(k, k + 1));

            // subtract original flux
            Conservative flux = Def::isHLLC ? HLLC(cells, face) : HLL(cells, face);
            cells.at(face.right).rezi -= cells.at(face.right).dt / cells.at(face.right).area * flux * face.length;

            // create new flux
            Conservative updatedFlux = Conservative(0, face.nx, face.ny, 0);
            double p_1 = Primitive::computePV(cells.at(k).w).p;
            double p_2 = Primitive::computePV(cells.at(k + Def::xCells).w).p;
            double p_w = 1.5 * p_1 - 0.5 * p_2;
            updatedFlux = updatedFlux * p_w;
            cells.at(face.right).rezi += cells.at(face.right).dt / cells.at(face.right).area * updatedFlux * face.length;
        }
    }
}

double Scheme::computeRezi(const std::unordered_map<int, Cell> &cells) {
    double res = 0;
    for (int i = 0; i < Def::inner; ++i) {
        int k = Def::innerIndex(i);
        res += pow(cells.at(k).rezi.r1 / cells.at(k).dt, 2) * cells.at(k).area;
        if (_isnan(res)) {
            cells.at(k).rezi.toString();
            Def::error = true;
            break;
        }
    }
    res = log(sqrt(res));
    return res;
}

void Scheme::updateCells(std::unordered_map<int, Cell> &cells) {
    for (int i = 0; i < Def::inner; ++i) {
        int k = Def::innerIndex(i);
        cells.at(k).w += cells.at(k).rezi;
        cells.at(k).rezi = 0;
    }
}

double Scheme::computeCP(double p_inner) {
    return (p_inner - Bound::p_infty) / (0.5 * Bound::rho_infty * pow(Bound::u_infty, 2) + pow(Bound::v_infty, 2));
}

Conservative Scheme::HLLC(const std::unordered_map<int, Cell> &cells, const Interface &face) {
    Conservative wl = cells.at(face.left).w;
    Conservative wr = cells.at(face.right).w;
    Primitive pvl = Primitive::computePV(wl);
    Primitive pvr = Primitive::computePV(wr);

    double ql = pvl.u * face.nx + pvl.v * face.ny; // normálová rychlost
    double qr = pvr.u * face.nx + pvr.v * face.ny; // DP - \tilde u

//    double q_bar = Scheme::bar(pvl.rho, pvr.rho, ql, qr); // viz Toro
    double h_bar = Scheme::bar(pvl.rho, pvr.rho, pvl.h, pvr.h);
    double u_bar = Scheme::bar(pvl.rho, pvr.rho, pvl.u, pvr.u);
    double v_bar = Scheme::bar(pvl.rho, pvr.rho, pvl.v, pvr.v);
    double U_bar_sq = pow(u_bar, 2) + pow(v_bar, 2); // certified J. Holman verze
    double c_bar = sqrt((Def::KAPPA - 1) * (h_bar - 0.5 * U_bar_sq));
    double q_bar = u_bar * face.nx + v_bar * face.ny;

    double lambda_1 = ql - pvl.c;
    double lambda_m = qr + pvr.c;
    double lambda_1Roe = q_bar - c_bar;
    double lambda_mRoe = q_bar + c_bar;

    double SL = fmin(lambda_1, lambda_1Roe);
    double SR = fmax(lambda_m, lambda_mRoe);
    double SM = (pvr.rho * qr * (SR - qr) - pvl.rho * ql * (SL - ql) + pvl.p - pvr.p) / (pvr.rho * (SR - qr) - pvl.rho * (SL - ql));

    double p_star = pvl.rho * (ql - SL) * (ql - SM) + pvl.p;

    Conservative wlStar = 1 / (SL - SM) * Scheme::fluxStar(face, wl, ql, SL, SM, pvl.p, p_star);
    Conservative wrStar = 1 / (SR - SM) * Scheme::fluxStar(face, wr, qr, SR, SM, pvr.p, p_star);

    if (SL > 0) {
        return Scheme::flux(face, wl, ql, pvl.p);
    } else if (SL <= 0 && 0 < SM) {
        return Scheme::flux(face, wlStar, SM, p_star);
    } else if (SM <= 0 && 0 <= SR) {
        return Scheme::flux(face, wrStar, SM, p_star);
    } else if (SR < 0) {
        return Scheme::flux(face, wr, qr, pvr.p);
    } else {
        std::cout << "Scheme::HLLC: unreachable state \n";
    }
}

Conservative Scheme::flux(Interface face, Conservative w, double q, double p) {
    Conservative res{};

    res.r1 = w.r1 * q;
    res.r2 = w.r2 * q + p * face.nx;
    res.r3 = w.r3 * q + p * face.ny;
    res.r4 = (w.r4 + p) * q;
    return res;
}

Conservative Scheme::fluxStar(Interface face, Conservative w, double q, double S, double SM, double p, double p_star) {
    Conservative res{};

    res.r1 = w.r1 * (S - q) + 0;
    res.r2 = w.r2 * (S - q) + (p_star - p) * face.nx;
    res.r3 = w.r3 * (S - q) + (p_star - p) * face.ny;
    res.r4 = w.r4 * (S - q) + p_star * SM - p * q;
    return res;
}

double Scheme::bar(double rho_l, double rho_r, double vl, double vr) {
    return (sqrt(rho_l) * vl + sqrt(rho_r) * vr) / (sqrt(rho_l) + sqrt(rho_r));
}