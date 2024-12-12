//
// Created by petrs on 01.12.2023.
//

#ifndef GAMM_PRIMITIVE_H
#define GAMM_PRIMITIVE_H


#include <iostream>
#include "Conservative.h"

class Primitive
{

public:
  double rho, rhoU, rhoE;
  double u, c, p, h;

  // Constructor
  Primitive () = default;

  Primitive (Conservative w);

  // Methods
  static Primitive computePV (Conservative w);
};


#endif //GAMM_PRIMITIVE_H
