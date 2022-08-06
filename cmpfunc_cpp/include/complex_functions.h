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

#include "./complex_operators.h"

long double inf_norm(const std::complex<long double> &z);
bool isfinite(const std::complex<long double> &z);

long double log_r(long double x);
long double dilog_r(long double x);

const std::complex<long double> II(0.0, 1.0);
std::complex<long double> log_c(std::complex<long double> z);
std::complex<long double> dilog_c(std::complex<long double> z);

std::complex<long double> log_c_angle(std::complex<long double> z,
                                      long double angle);
std::complex<long double> expm1(const std::complex<long double> &z);
std::complex<long double> log1p(const std::complex<long double> &z);

std::complex<long double> LogGamma(std::complex<long double> z);
std::complex<long double> Gamma_inv(const std::complex<long double> &z);

std::complex<long double> LBesselJ(long double k, std::complex<long double> z);
std::complex<long double> LBesselK(long double k, std::complex<long double> z);
std::complex<long double> CBesselK(std::complex<long double> nu,
                                   std::complex<long double> z,
                                   bool warn = true);

std::complex<long double> Hyp2F1(std::complex<long double> a,
                                 std::complex<long double> b,
                                 std::complex<long double> c,
                                 std::complex<long double> z);

std::complex<long double> Hyp2F1(std::complex<long double> a,
                                 std::complex<long double> b,
                                 std::complex<long double> c,
                                 std::complex<long double> z);

std::complex<long double> incBeta(std::complex<long double> x,
                                  std::complex<long double> a,
                                  std::complex<long double> b);

std::complex<long double> incGamma2(std::complex<long double> x,
                                    std::complex<long double> a);

std::complex<long double> incGamma3(std::complex<long double> x,
                                    std::complex<long double> a,
                                    std::complex<long double> b);
