//
// Created by petrs on 01.12.2023.
//

#include <cmath>
#include "Primitive.h"

constexpr double KAPPA = 1.4;

Primitive Primitive::computePV(Conservative w) {
    Primitive res{};

    res.rho = w.r1;
    res.rhoU = w.r2;
    res.rhoE = w.r3;

    res.u = w.r2 / w.r1;
    res.p = (KAPPA - 1) * (res.rhoE - 0.5 * res.rho * pow(res.u, 2));

    res.c = sqrt((KAPPA*res.p) / res.rho);
    res.h = (KAPPA * res.p) / (res.rho * (KAPPA - 1));

    return res;
} //kontrolováno - přestos tam měl chybu blbče
