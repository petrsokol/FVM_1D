//
// Created by petrs on 15.10.2024.
//

#include <cstdlib>
#include <vector>
#include <cmath>
#include <fstream>
#include "structures/Conservative.h"
#include "structures/Primitive.h"

/**
 * solution to problem no. 6 from
 * https://moodle-vyuka.cvut.cz/course/view.php?id=12877
 */


const std::string path = R"(C:\Users\petrs\Documents\CTU\BP\Charts\Data)";

// numerical scheme computation
/*------------------------------------------------------------*/
double computeTimeStep (double CFL, const std::vector<Conservative> &w, double *t, double dx);

double computeRezi (const std::vector<Conservative> &w, const std::vector<Conservative> &wn, double dt, double dx);

double bar (double rho_l, double rho_r, double vl, double vr);

// HLL scheme + support functions
Conservative HLL (Conservative wl, Conservative wr);

Conservative flux (Conservative w, double q, double p);

// HLLC scheme + support functions
Conservative HLLC (Conservative wl, Conservative wr);

Conservative fluxStar (Conservative w, double q, double S, double SM, double p, double p_star);

// second order limiter
Conservative minmod (Conservative a, Conservative b);

double minmod (double a, double b);

/*------------------------------------------------------------*/

void exportData (const std::string &path, const std::string &filename, const std::vector<Conservative> &w, double dx);

std::string fileName (bool is_hll, bool second_order, int snapshotCount, int innerCells);

// implementation
/*------------------------------------------------------------*/
double setInitialConditions (std::vector<Conservative> &w, std::vector<Conservative> &wn, Conservative W_L,
                             Conservative W_R, double midBound);

void updateBounds (const std::vector<Conservative> &w, std::vector<Conservative> &wn);

void computeScheme (const std::vector<Conservative> &w, std::vector<Conservative> &wn,
                    double dt, bool isHLL, bool isSecOrd, double dx);

void runExperiment (int innerCells, bool isHLL, bool isSecOrd, double T);

/*------------------------------------------------------------*/

int getGhostLayers ();

double getKappa ();

Conservative setWithRhoUP (double rho, double u, double p);

int main ()
{
  int innerCells = 500;
  double T = 0.2;
  runExperiment(innerCells, true, false,  T);
  runExperiment(innerCells, false, false, T);
  runExperiment(innerCells, true, true,   T);
  runExperiment(innerCells, false, true,  T);

  // Define the runPythonScript to execute the Python script
  const char *runPythonScript = R"(python "C:\Users\petrs\Documents\CTU\BP\PYTHON-scripts\animateChart.py")";

  // Check if the system runPythonScript succeeded
  if (system(runPythonScript) == -1) {
    printf("Failed to run Python script.\n");
    return 1;
  } else {
    printf("Python script executed successfully.\n");
  }

  return EXIT_SUCCESS;
  /*
   * todo
   *   - automatickej output tak, abych mohl sestavit grafy v různých časových intervalech
   *   - minmod a tvd rekonstrukce druhýho řádu přesnosti
   *   - přepínání na libovolnou konfiguraci schématu a řádu přesnosti
   */
}

void runExperiment (const int innerCells, const bool isHLL, const bool isSecOrd, const double T)
{
  // declare vectors for two consecutive time steps
  std::vector<Conservative> w(innerCells + 2 * getGhostLayers());
  std::vector<Conservative> wn(innerCells + 2 * getGhostLayers());

  // set initial condition
  Conservative W_L = setWithRhoUP(1, 0, 1);
  Conservative W_R = setWithRhoUP(0.1, 0, 0.1795);
  double dx = setInitialConditions(w, wn, W_L, W_R, 1);

  // start the process
  double t = 0;
  while (t < T) {
    double dt = computeTimeStep(0.4, w, &t, dx);

    computeScheme(w, wn, dt, isHLL, isSecOrd, dx);

    updateBounds(w, wn);

    // update the flux-conservative values for each cell
    w = wn;
  }

  // export data to a external file
  exportData(path, fileName(isHLL, isSecOrd, 20, innerCells), w, dx);
}

Conservative setWithRhoUP (const double rho, const double u, const double p)
{
  const double KAPPA = getKappa();
  double rhoU = rho * u;
  double rhoE = p / (KAPPA - 1) + 0.5 * rho * u * u;
  return Conservative(rho, rhoU, rhoE);
}

void computeScheme (const std::vector<Conservative> &w, std::vector<Conservative> &wn,
                    const double dt, const bool isHLL, const bool isSecOrd, const double dx)
{
  const int GHOST_LAYERS = getGhostLayers();
  const int firstInner = GHOST_LAYERS;
  const int innerCells = (int) w.size() - 2 * GHOST_LAYERS;
  for (int i = firstInner - 1; i < firstInner + innerCells; ++i) {
    int l = i;
    int r = i + 1;
    int so = isSecOrd ? 1 : 0;

    Conservative sigma_l_dopr = (w.at(r) - w.at(l)) / dx;
    Conservative sigma_l_zpet = (w.at(l) - w.at(l - 1)) / dx;
    Conservative sigma_r_dopr = (w.at(r + 1) - w.at(r)) / dx;
    Conservative sigma_r_zpet = (w.at(r) - w.at(l)) / dx; // = sigma_l_dopr

    Conservative sigma_l = minmod(sigma_l_dopr, sigma_l_zpet) * so;
    Conservative sigma_r = minmod(sigma_r_dopr, sigma_r_zpet) * so;

    Conservative wl = w.at(l) + dx / 2 * sigma_l;
    Conservative wr = w.at(r) - dx / 2 * sigma_r;

    Conservative flux = isHLL ? HLL(wl, wr) : HLLC(wl, wr);
    wn.at(l) -= dt / dx * flux;
    wn.at(r) += dt / dx * flux;
  }
}

double minmod (double a, double b)
{
  if (a * b <= 0) {
    return 0;
  } else if (fabs(a) <= fabs(b) && a * b > 0) {
    return a;
  } else if (fabs(a) > fabs(b) && a * b > 0) {
    return b;
  } else {
    printf("sus\n");
    return 0;
  }
}

Conservative minmod (Conservative a, Conservative b)
{
  Conservative res;
  res.r1 = minmod(a.r1, b.r1);
  res.r2 = minmod(a.r2, b.r2);
  res.r3 = minmod(a.r3, b.r3);
  return res;
}

void updateBounds (const std::vector<Conservative> &w, std::vector<Conservative> &wn)
{
  const int cells = (int) w.size();
  const int limit = getGhostLayers();
  for (int i = 0; i < limit; ++i) {
    wn.at(i) = w.at(limit);
    wn.at(cells - 1 - i) = w.at(cells - 1 - limit);
  }
}

double setInitialConditions (std::vector<Conservative> &w, std::vector<Conservative> &wn,
                             const Conservative W_L, const Conservative W_R, const double midBound)
{
  // returns value of dx;
  constexpr double X_LOWER_BOUND = 0;
  constexpr double X_UPPER_BOUND = 2;

  const int GHOST_LAYERS = getGhostLayers();
  const int innerCells = (int) w.size() - 2 * GHOST_LAYERS;
  double dx = (X_UPPER_BOUND - X_LOWER_BOUND) / innerCells;

  int cells = (int) w.size();

  for (int i = 0; i < cells; ++i) {
    double x = (i - GHOST_LAYERS - 1) * dx;
    if (x < midBound) {
      w.at(i) = W_L;
    } else {
      w.at(i) = W_R;
    }
  }

  wn = w;

  return dx;
}


std::string fileName (bool is_hll, bool second_order, int snapshotCount, int innerCells)
{
  std::string scheme = is_hll ? "HLL" : "HLLC";
  std::string order = second_order == 1 ? "higherOrder" : "firstOrder";
  std::string slice = std::to_string(snapshotCount);
  std::string innerCellCount = std::to_string(innerCells);
  std::string name = order + "_" + scheme + "_" + slice + "_" + innerCellCount;
  return name;
}

double computeTimeStep (double CFL, const std::vector<Conservative> &w, double *t, const double dx)
{
  // problematic if all initial velocities are zero
  double res = 1e-4;
  const int cells = (int) w.size();
  for (int i = 0; i < cells; ++i) {
    auto pv = Primitive(w.at(i));
    res = fmin(res, dx / (fabs(pv.u * pv.c)));
  }
  res = res * CFL;
  *t += res;
//  printf("computeTimeStep: %f\n", res);
  return res;
}

double computeRezi (const std::vector<Conservative> &w, const std::vector<Conservative> &wn, double dt, const double dx)
{
  double res = 0;
  const int cells = (int) w.size();
  for (int i = 0; i < cells; ++i) {
    res += pow((wn.at(i).r1 - w.at(i).r1) / dt, 2) * dx;
  }
  return log10(sqrt(res));
}

Conservative HLLC (Conservative wl, Conservative wr)
{
  Conservative res;
  auto pvl = Primitive(wl);
  auto pvr = Primitive(wr);
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
  double SM = (pvr.rho * qr * (SR - qr) - pvl.rho * ql * (SL - ql) + pvl.p - pvr.p) /
              (pvr.rho * (SR - qr) - pvl.rho * (SL - ql));

  double p_star = pvl.rho * (ql - SL) * (ql - SM) + pvl.p;

  Conservative wlStar = 1 / (SL - SM) * fluxStar(wl, ql, SL, SM, pvl.p, p_star);
  Conservative wrStar = 1 / (SR - SM) * fluxStar(wr, qr, SR, SM, pvr.p, p_star);

  if (SL > 0) {
    res = flux(wl, ql, pvl.p);
  } else if (SL <= 0 && 0 < SM) {
    res = flux(wlStar, SM, p_star);
  } else if (SM <= 0 && 0 <= SR) {
    res = flux(wrStar, SM, p_star);
  } else if (SR < 0) {
    res = flux(wr, qr, pvr.p);
  } else {
    std::cout << "Scheme::HLLC: unreachable state \n";
  }

  return res;
}

double bar (double rho_l, double rho_r, double vl, double vr)
{
  return (sqrt(rho_l) * vl + sqrt(rho_r) * vr) / (sqrt(rho_l) + sqrt(rho_r));
}

Conservative HLL (Conservative wl, Conservative wr)
{
  Conservative res{};

  /* compute primitive variables */
  auto pvl = Primitive(wl);
  auto pvr = Primitive(wr);

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

Conservative fluxStar (Conservative w, double q, double S, double SM, double p, double p_star)
{
  Conservative res{};

  res.r1 = w.r1 * (S - q) + 0;
  res.r2 = w.r2 * (S - q) + (p_star - p);
  res.r3 = w.r3 * (S - q) + p_star * SM - p * q;

  return res;
}

Conservative flux (Conservative w, double q, double p)
{
  Conservative res{};

  res.r1 = w.r1 * q;
  res.r2 = w.r2 * q + p;
  res.r3 = (w.r3 + p) * q;

  return res;
}

void exportData (const std::string &filePath, const std::string &fileName, const std::vector<Conservative> &w,
                 const double dx)
{
  std::ofstream output(filePath + "\\" + fileName + ".dat");
  const int innerCells = (int) w.size() - 2 * getGhostLayers();
  for (int i = 0; i < innerCells; ++i) {
    output << i * dx << " " << w[i].r1 << " " << w[i].r2 << " " << w[i].r3 << std::endl;
  }
  output.close();
  std::cout << "Exported successfully" << std::endl;
}

int getGhostLayers ()
{
  return 2;
}

double getKappa ()
{
  return 1.4;
}