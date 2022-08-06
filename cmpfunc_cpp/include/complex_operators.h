/*
 * =====================================================================================
 *
 *       Filename:  complex_functions.h
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

#pragma once

#include <complex>
#include <complex_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_zeta.h>
#include <iostream>
#include <numeric>
#include <vector>

using namespace sp_bessel;

std::complex<long double> operator*(const std::complex<long double> &a,
                                    const double &b);
std::complex<long double> operator*(const double &a,
                                    const std::complex<long double> &b);
std::complex<long double> operator/(const std::complex<long double> &a,
                                    const double &b);
std::complex<long double> operator/(const double &a,
                                    const std::complex<long double> &b);
std::complex<long double> operator-(const std::complex<long double> &a,
                                    const double &b);
std::complex<long double> operator-(const double &a,
                                    const std::complex<long double> &b);
std::complex<long double> operator+(const std::complex<long double> &a,
                                    const double &b);
std::complex<long double> operator+(const double &a,
                                    const std::complex<long double> &b);

bool operator==(const std::complex<long double> &z, const double &a);
bool operator!=(const std::complex<long double> &z, const double &a);
bool operator==(const double &a, const std::complex<long double> &z);
bool operator!=(const double &a, const std::complex<long double> &z);

std::complex<long double> operator+(const std::complex<long double> &z,
                                    const int n);
std::complex<long double> operator-(const std::complex<long double> &z,
                                    const int n);
std::complex<long double> operator*(const std::complex<long double> &z,
                                    const int n);
std::complex<long double> operator/(const std::complex<long double> &z,
                                    const int n);
std::complex<long double> operator+(const int n,
                                    const std::complex<long double> &z);
std::complex<long double> operator-(const int n,
                                    const std::complex<long double> &z);
std::complex<long double> operator*(const int n,
                                    const std::complex<long double> &z);
std::complex<long double> operator/(const int n,
                                    const std::complex<long double> &z);
std::complex<long double> operator+(const std::complex<long double> &z,
                                    const unsigned int n);
std::complex<long double> operator-(const std::complex<long double> &z,
                                    const unsigned int n);
std::complex<long double> operator*(const std::complex<long double> &z,
                                    const unsigned int n);
std::complex<long double> operator/(const std::complex<long double> &z,
                                    const unsigned int n);
std::complex<long double> operator+(const unsigned int n,
                                    const std::complex<long double> &z);
std::complex<long double> operator-(const unsigned int n,
                                    const std::complex<long double> &z);
std::complex<long double> operator*(const unsigned int n,
                                    const std::complex<long double> &z);
std::complex<long double> operator/(const unsigned int n,
                                    const std::complex<long double> &z);
std::complex<long double>
operator*(const std::vector<std::complex<long double>> &c1,
          const std::vector<std::complex<long double>> &c2);

bool operator==(const std::complex<long double> &z, const int n);
bool operator!=(const std::complex<long double> &z, const int n);
bool operator==(const int n, const std::complex<long double> &z);
bool operator!=(const int n, const std::complex<long double> &z);
bool operator==(const std::complex<long double> &z, const unsigned int n);
bool operator!=(const std::complex<long double> &z, const unsigned int n);
bool operator==(const unsigned int n, const std::complex<long double> &z);
bool operator!=(const unsigned int n, const std::complex<long double> &z);

std::complex<long double> pow(const long double &a,
                              const std::complex<long double> &b);
std::complex<long double> pow(const std::complex<long double> &a,
                              const long double &b);
std::complex<long double> pow(const std::complex<long double> &a,
                              const std::complex<long double> &b);
