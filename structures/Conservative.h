//
// Created by petrs on 01.12.2023.
//

#ifndef GAMM_CONSERVATIVE_H
#define GAMM_CONSERVATIVE_H


#include <iostream>

class Conservative
{

public:
  // parameters
  double r1, r2, r3;

  // Constructor
  explicit Conservative (double r1 = 0.0, double r2 = 0.0, double r3 = 0.0) : r1(r1), r2(r2), r3(r3) {}

  // Methods
  void toString () const
  {
    std::cout << "[" << r1 << ";" << r2 << ";" << r3 << "]\n";
  }

  // Overload the + operator for vector addition
  Conservative operator+ (const Conservative &other) const
  {
    return Conservative(r1 + other.r1, r2 + other.r2, r3 + other.r3);
  }

  // Overload the - operator for vector subtraction
  Conservative operator- (const Conservative &other) const
  {
    return Conservative(r1 - other.r1, r2 - other.r2, r3 - other.r3);
  }

  // Overload the = operator for whatever
  Conservative operator= (double scalar)
  {
    return Conservative(r1 = scalar, r2 = scalar, r3 = scalar);
  }

  // Overload the += operator
  Conservative &operator+= (const Conservative &other)
  {
    // Perform the addition for each member
    r1 += other.r1;
    r2 += other.r2;
    r3 += other.r3;

    // Return the modified object
    return *this;
  }

  // Overload the -= operator
  Conservative &operator-= (const Conservative &other)
  {
    // Perform the addition for each member
    r1 -= other.r1;
    r2 -= other.r2;
    r3 -= other.r3;

    // Return the modified object
    return *this;
  }

  // Overload the * operator for scalar multiplication
  Conservative operator* (double scalar) const
  {
    return Conservative(r1 * scalar, r2 * scalar, r3 * scalar);
  }

  // Overload the / operator for scalar division
  Conservative operator/ (double scalar) const
  {
    // Check for division by zero
    if (scalar == 0.0) {
      std::cerr << "Error: Division by zero\n";
      // You might want to handle this error differently based on your needs
      exit(EXIT_FAILURE);
    }
    return Conservative(r1 / scalar, r2 / scalar, r3 / scalar);
  }
};

Conservative operator* (double scalar, const Conservative &vec);


#endif //GAMM_CONSERVATIVE_H
