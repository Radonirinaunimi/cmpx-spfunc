/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Test library that computes special functions wtih cmpx args
 *
 *        Version:  1.0
 *        Created:  08/06/2022 23:16:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "../cmpfunc_cpp/include/complex_functions.h"

struct err_message : public std::exception {
  const char *what() const throw() { return "Wrong Parameters!!"; }
};

int main() {

  std::complex<long double> a(1, -4);
  std::complex<long double> b(-2, 2);
  std::complex<long double> c(3, -9);
  std::complex<long double> z(0, 26);

  std::cout << std::scientific;
  // Compute the Hypergeometric 2F1 with complex arguments
  std::complex<long double> hyp2f1 = Hyp2F1(a, b, c, z);

  std::cout << "Hyp2F1: "
            << "real=" << hyp2f1.real() << "; img=" << hyp2f1.imag()
            << std::endl;

  // Compute the incomplete Beta with complex arguments
  std::complex<long double> inbeta = incBeta(a, b, c);

  std::cout << "InBeta: "
            << "real=" << inbeta.real() << "; img=" << inbeta.imag()
            << std::endl;

  return 0;
}
