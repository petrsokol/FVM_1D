//
// Created by petrs on 01.12.2023.
//

#include <cmath>
#include "Primitive.h"

constexpr double KAPPA = 1.4;

// Constructor definition
Primitive::Primitive(Conservative w)
{
  rho = w.r1;
  rhoU = w.r2;
  rhoE = w.r3;
  u = w.r2 / w.r1;
  p = (KAPPA - 1) * (rhoE - 0.5 * rho * pow(u, 2));
  c = sqrt((KAPPA*p) / rho);
  h = (KAPPA * p) / (rho * (KAPPA - 1));
}
