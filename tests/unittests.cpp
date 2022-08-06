/*
 * =====================================================================================
 *
 *       Filename:  TestDY.cpp
 *
 *    Description:  This file performs the benchmark of the full hadronic
 * expression and all the expressions it depends on with the HqT code
 *                  http://theory.fi.infn.it/grazzini/codes.html.
 *
 *        Version:  1.0
 *        Created:  18/02/2021 00:15:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#define CATCH_CONFIG_MAIN

#include "../cmpfunc_cpp/include/complex_functions.h"
#include "catch.hpp"

// Define Global Variables
std::complex<long double> A(1, -4);
std::complex<long double> B(-2, 2);
std::complex<long double> C(3, -9);
std::complex<long double> Z(0, 26);

TEST_CASE("Test Hypergeometric function 2F1") {
  std::complex<long double> hyp2f1 = Hyp2F1(A, B, C, Z);

  long double real_value = hyp2f1.real();
  long double imag_value = hyp2f1.imag();

  long double real_refer = -2.427077e+00;
  long double imag_refer = -6.856067e+00;

  REQUIRE(real_value == Approx(real_refer).epsilon(1e-5));
  REQUIRE(imag_value == Approx(imag_refer).epsilon(1e-5));
}
