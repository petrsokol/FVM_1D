#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "structures/Conservative.h"
#include "structures/Primitive.h"

using namespace std;

const double KAPPA = 1.4;
double xLowerBound, xUpperBound, dx;
double t, T, dt;
int innerCells;
//
//struct conservative {
//    double r1, r2, r3;
//}; //vektory jako W, F
//
//struct primitive {
//    double rho, rhoU, rhoE, u, p, c, h;
//}; //pomocný vektor s fys. veličinami

void setTimeInterval(double lowerBound, double upperBound) {
    t = lowerBound;
    T = upperBound;
}

void setSpaceInterval(double lowerBound, double upperBound, int n) {
    xLowerBound = lowerBound;
    xUpperBound = upperBound;
    innerCells = n;
    innerCells = innerCells + 2;
    dx = (xUpperBound - xLowerBound) / innerCells;
}

void setInitialConditions(vector<double>& v, vector<Conservative>& values, Conservative vL, Conservative vR, double boundary) {
    for (int i = 0; i < innerCells; ++i) {
        v[i] = dx*(i-1);
    }
    for (int i = 0; i < innerCells; ++i) {
        if (v[i] <= boundary) { //pro datový typ double je třeba (v[i] <= boundary+1e-14)
            values[i] = vL;
        } else {
            values[i] = vR;
        }
    }
}

void exportData(const string& path, const string& filename, vector<Conservative> V, vector<double> x) {
    ofstream output(path + "\\" + filename + ".dat");
    for (int i = 1; i < innerCells; ++i) {
        output << x[i] << " " << V[i].r1 << endl;
    }
    output.close();
    cout << "Exported successfully" << endl;
}

void updateBoundary(vector<Conservative>& W) {
    W[0] = W[1];
    W[innerCells - 1] = W[innerCells - 2];
}
double computeTimeStep(double CFL, vector<Primitive> vars) {
    double timeStep = dx / (abs(vars[1].u) * vars[1].c);
    for (int i = 2; i < innerCells+1; ++i) {
        timeStep = min(timeStep, dx / (abs(vars[i].u) * vars[i].c));
    }
    return timeStep * CFL;
}

vector<Primitive> computePhysicalVariables(vector<Conservative> W) {
    vector<Primitive> vec(innerCells);
    double u,p;
    for (int i = 0; i < innerCells; ++i) {
        u = W[i].r2/W[i].r1;
        p = (KAPPA - 1)*(W[i].r3 - 0.5*W[i].r2*u);

        vec[i].rho = W[i].r1;
        vec[i].rhoU = W[i].r2;
        vec[i].rhoE = W[i].r3;
        vec[i].u = u;
        vec[i].p = p;
        vec[i].c = sqrt((KAPPA*p) / W[i].r1);
        vec[i].h = u + p / W[i].r1;
    }
    return vec;
}

vector<Conservative> computeFlux(vector<Primitive> pv) {
    vector<Conservative> f(innerCells);
    for (int i = 0; i < innerCells; ++i) {
        f[i].r1 = pv[i].rhoU;
        f[i].r2 = pv[i].rhoU * pv[i].u + pv[i].p;
        f[i].r3 = (pv[i].rhoE + pv[i].p) * pv[i].u;
    }
    return f;
}

Conservative fluxHLL(vector<Conservative> f, vector<Conservative> w, vector<Primitive> pv, int i) {
    Conservative res{};
    double SL = min(pv[i].u - pv[i].c, pv[i+1].u - pv[i+1].c);
    double SR = max(pv[i].u + pv[i].c, pv[i+1].u + pv[i+1].c);

    if (SL > 0) {
        res = f[i];
    } else if (SR < 0) {
        res = f[i+1];
    } else {
        res.r1 = (SR*f[i].r1 - SL*f[i+1].r1 + SL*SR*(w[i+1].r1 - w[i].r1))/(SR - SL);
        res.r2 = (SR*f[i].r2 - SL*f[i+1].r2 + SL*SR*(w[i+1].r2 - w[i].r2))/(SR - SL);
        res.r3 = (SR*f[i].r3 - SL*f[i+1].r3 + SL*SR*(w[i+1].r3 - w[i].r3))/(SR - SL);
    }
    return res;
}

double bar(double varL, double varR, vector<Primitive> pv, int i) {
    return (sqrt(pv[i].rho) * varL + sqrt(pv[i+1].rho) * varR) / (sqrt(pv[i].rho) + sqrt(pv[i+1].rho));
}

Conservative star(Conservative w, double u, double p_star, double p, double S, double SM) {
    Conservative w_star{};
    double omega = 1/(S-SM);
    w_star.r1 = omega * (w.r1*(S-u));
    w_star.r2 = omega * (w.r2*(S-u) + (p_star - p)); //p_star - p //absolutní hodnota spraví?
    w_star.r3 = omega * (w.r3*(S-u) - p*u + p_star*SM);
    return w_star;
}

Conservative fluxHLLC(vector<Conservative> f, vector<Conservative> w, vector<Primitive> pv, int i) {
    Conservative res{};
    int l = i;
    int r = i+1;
    double u_bar, c_bar, h_bar, q_barSq;
    double lambda_1, lambda_1Roe, lambda_m, lambda_mRoe;
    double SL, SR, SM;
    double p_star;
    Conservative w_star{};

    //TODO pomocí pointerů:
    double uL = pv[l].u;
    double cL = pv[l].c;
    double pL = pv[l].p;
    double hL = pv[l].h;
    double rhoL = pv[l].rho;

    double uR = pv[r].u;
    double cR = pv[r].c;
    double pR = pv[r].p;
    double hR = pv[r].h;
    double rhoR = pv[r].rho;

    u_bar = bar(uL, uR, pv, i);
    h_bar = bar(hL, hR, pv, i);
    q_barSq = bar(pow(uL, 2), pow(uR, 2), pv, i);
    c_bar = sqrt((KAPPA - 1)*(h_bar - q_barSq/2));

    lambda_1 = uL - cL;
    lambda_m = uR + cR;
    lambda_1Roe = u_bar - c_bar;
    lambda_mRoe = u_bar + c_bar;

    SL = min(lambda_1, lambda_1Roe);
    SR = max(lambda_m, lambda_mRoe);
    SM = (rhoR*uR*(SR-uR) - rhoL*uL*(SL-uL) + pL - pR) /
            (rhoR*(SR-uR) - rhoL*(SL-uL));
    p_star = rhoL*(uL - SL)*(uL - SM) + pL; //TODO jaktože platí rovnost (2.97)?

    if (SL > 0) {
        res = f[l];
    } else if (SL <= 0 && 0 < SM) {
        w_star = star(w[l], uL, p_star, pL, SL, SM);
        res.r1 = SM * (w_star.r1);
        res.r2 = SM * (w_star.r2) + p_star;
        res.r3 = SM * (w_star.r3 + p_star);

    } else if (SM <= 0 && 0 <= SR) {
        w_star = star(w[r], uR, p_star, pR, SR, SM);
        res.r1 = SM * (w_star.r1); //TODO opakuješ se
        res.r2 = SM * (w_star.r2) + p_star;
        res.r3 = SM * (w_star.r3 + p_star);
    } else {
        res = f[r]; // todo what if NaN or inf??
    }
    return res;
}

int not_main() {
    string path = R"(C:\Users\petrs\Documents\CTU\BP\Charts\Data)";
    setSpaceInterval(0, 1.5, 150);
    setTimeInterval(0, 0.2);

    vector<Conservative> w(innerCells);
    vector<Conservative> wN(innerCells);
    vector<Conservative> f(innerCells);

    Conservative flux{};

    vector<Primitive> PV(innerCells);
    vector<double> x(innerCells);

    Conservative wL = Conservative(0.445, 0.311, 8.9280);
    Conservative wR = Conservative(0.5, 0, 1.4275);
    setInitialConditions(x, w, wL, wR, 0.7);

    while(t < T) {
        PV = computePhysicalVariables(w);
        dt = computeTimeStep(0.8, PV);
        t += dt;
        printf("dt = %f\n", dt);
        f = computeFlux(PV);
        wN = w;

        for (int i = 0; i < innerCells - 1; ++i) {
            flux = fluxHLL(f, w, PV, i); //buď fluxHLL(...) nebo fluxHLLC(...)
            wN[i] -= dt / dx * flux;
            wN[i+1] += dt / dx * flux;
        }

        updateBoundary(wN);
        w = wN;
    }
    exportData(path, "hll_1D_ref_2", w, x);
    cout << "final time: " << T << endl;
}