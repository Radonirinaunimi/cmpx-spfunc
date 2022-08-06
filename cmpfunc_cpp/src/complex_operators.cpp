/*
 * =====================================================================================
 *
 *       Filename:  complex_functions.cpp
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

#include "../include/complex_operators.h"

/*#####################################################################################
 *             DEFINITION OF COMPLEX OPERATORS AND FUNCTIONS
 *####################################################################################*/

std::complex<long double> operator*(const std::complex<long double> &a,
                                    const double &b) {
  return a * static_cast<long double>(b);
}

std::complex<long double> operator*(const double &a,
                                    const std::complex<long double> &b) {
  return static_cast<long double>(a) * b;
}

std::complex<long double> operator/(const std::complex<long double> &a,
                                    const double &b) {
  return a / static_cast<long double>(b);
}

std::complex<long double> operator/(const double &a,
                                    const std::complex<long double> &b) {
  return static_cast<long double>(a) / b;
}

std::complex<long double> operator-(const std::complex<long double> &a,
                                    const double &b) {
  return a - static_cast<long double>(b);
}

std::complex<long double> operator-(const double &a,
                                    const std::complex<long double> &b) {
  return static_cast<long double>(a) - b;
}

std::complex<long double> operator+(const std::complex<long double> &a,
                                    const double &b) {
  return a + static_cast<long double>(b);
}

std::complex<long double> operator+(const double &a,
                                    const std::complex<long double> &b) {
  return static_cast<long double>(a) + b;
}

bool operator==(const std::complex<long double> &z, const double &a) {
  return (z == static_cast<long double>(a));
}
bool operator!=(const std::complex<long double> &z, const double &a) {
  return (z != static_cast<long double>(a));
}
bool operator==(const double &a, const std::complex<long double> &z) {
  return (static_cast<long double>(a) == z);
}
bool operator!=(const double &a, const std::complex<long double> &z) {
  return (static_cast<long double>(a) != z);
}

std::complex<long double> operator+(const std::complex<long double> &z,
                                    const int n) {
  return (z + static_cast<long double>(n));
}

std::complex<long double> operator-(const std::complex<long double> &z,
                                    const int n) {
  return (z - static_cast<long double>(n));
}

std::complex<long double> operator*(const std::complex<long double> &z,
                                    const int n) {
  return (z * static_cast<long double>(n));
}

std::complex<long double> operator/(const std::complex<long double> &z,
                                    const int n) {
  return (z / static_cast<long double>(n));
}

std::complex<long double> operator+(const int n,
                                    const std::complex<long double> &z) {
  return (static_cast<long double>(n) + z);
}

std::complex<long double> operator-(const int n,
                                    const std::complex<long double> &z) {
  return (static_cast<long double>(n) - z);
}

std::complex<long double> operator*(const int n,
                                    const std::complex<long double> &z) {
  return (static_cast<long double>(n) * z);
}

std::complex<long double> operator/(const int n,
                                    const std::complex<long double> &z) {
  return (static_cast<long double>(n) / z);
}

// Overload Operators
std::complex<long double> operator+(const std::complex<long double> &z,
                                    const unsigned int n) {
  return (z + static_cast<long double>(n));
}

std::complex<long double> operator-(const std::complex<long double> &z,
                                    const unsigned int n) {
  return (z - static_cast<long double>(n));
}

std::complex<long double> operator*(const std::complex<long double> &z,
                                    const unsigned int n) {
  return (z * static_cast<long double>(n));
}

std::complex<long double> operator/(const std::complex<long double> &z,
                                    const unsigned int n) {
  return (z / static_cast<long double>(n));
}

std::complex<long double> operator+(const unsigned int n,
                                    const std::complex<long double> &z) {
  return (static_cast<long double>(n) + z);
}

std::complex<long double> operator-(const unsigned int n,
                                    const std::complex<long double> &z) {
  return (static_cast<long double>(n) - z);
}

std::complex<long double> operator*(const unsigned int n,
                                    const std::complex<long double> &z) {
  return (static_cast<long double>(n) * z);
}

std::complex<long double> operator/(const unsigned int n,
                                    const std::complex<long double> &z) {
  return (static_cast<long double>(n) / z);
}

bool operator==(const std::complex<long double> &z, const int n) {
  return (z == static_cast<long double>(n));
}

bool operator!=(const std::complex<long double> &z, const int n) {
  return (z != static_cast<long double>(n));
}

bool operator==(const int n, const std::complex<long double> &z) {
  return (static_cast<long double>(n) == z);
}

bool operator!=(const int n, const std::complex<long double> &z) {
  return (static_cast<long double>(n) != z);
}

bool operator==(const std::complex<long double> &z, const unsigned int n) {
  return (z == static_cast<long double>(n));
}

bool operator!=(const std::complex<long double> &z, const unsigned int n) {
  return (z != static_cast<long double>(n));
}

bool operator==(const unsigned int n, const std::complex<long double> &z) {
  return (static_cast<long double>(n) == z);
}

bool operator!=(const unsigned int n, const std::complex<long double> &z) {
  return (static_cast<long double>(n) != z);
}

std::complex<long double> pow(const long double &a,
                              const std::complex<long double> &b) {
  std::complex<long double> A(a, 0.);
  return (std::exp(b * std::log(A)));
}

std::complex<long double> pow(const std::complex<long double> &a,
                              const long double &b) {
  std::complex<long double> B(b, 0.);
  return (std::exp(B * std::log(a)));
}
std::complex<long double> pow(const std::complex<long double> &a,
                              const std::complex<long double> &b) {
  return (std::exp(b * std::log(a)));
}

long double inf_norm(const std::complex<long double> &z) {
  return std::max(std::abs(std::real(z)), std::abs(std::imag(z)));
}

bool isfinite(const std::complex<long double> &z) {
  const long double x = std::real(z), y = std::imag(z);
  return (std::isfinite(x) && std::isfinite(y));
}

// Inner Product
std::complex<long double>
operator*(const std::vector<std::complex<long double>> &c1,
          const std::vector<std::complex<long double>> &c2) {
  std::complex<long double> result;
  std::complex<long double> zero(0., 0.);
  result = std::inner_product(c1.begin(), c1.end(), c2.begin(), zero);

  return result;
}
