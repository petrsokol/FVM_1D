//
// Created by petrs on 01.12.2023.
//

#include "Conservative.h"

Conservative operator*(double scalar, const Conservative& vec) {
    return vec * scalar;
}