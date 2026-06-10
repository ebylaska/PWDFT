#ifndef _UNITS_HPP_
#define _UNITS_HPP_

#pragma once // Modern, clean, and avoids redundancy

#include <numbers> // C++20: The gold standard for math constants

namespace pwdft::units {

// ----------------------------
// Fundamental constants (SI)
// ----------------------------
// 'inline' ensures one single instance across the whole program
inline constexpr double KB = 1.380649e-23;          // J/K
inline constexpr double H  = 6.62607015e-34;        // J·s
inline constexpr double PI = 3.14159265358979323846; 
inline constexpr double C    = 2.99792458e10;       // cm/s
inline constexpr double C_SI = 2.99792458e8;        // m/s
inline constexpr double C_CM = 2.99792458e10;       // cm/s

// ----------------------------
// Thermodynamic constants
// ----------------------------
inline constexpr double R = 8.31446261815324e-3;    // kJ/mol·K

// ----------------------------
// Spectroscopic conversions
// ----------------------------
// These are computed at compile-time via "constant folding"
inline constexpr double HC_OVER_K = (H * C) / KB;   // cm·K
inline constexpr double AU_FREQ_TO_CM = 219474.6313705;

// ----------------------------
// Mass conversions
// ----------------------------
inline constexpr double AMU_TO_ME = 1822.888486209;
inline constexpr double ME_TO_AMU = 1.0 / AMU_TO_ME;
inline constexpr double AMU_TO_KG = 1.66053906660e-27;
inline constexpr double ME_TO_KG = 9.1093837015e-31;


// ----------------------------
// Length conversions
// ----------------------------
inline constexpr double ANG2_TO_M2 = 1.0e-20;
inline constexpr double BOHR_TO_M = 5.29177210903e-11;
inline constexpr double BOHR_TO_ANG = 0.529177210903;
inline constexpr double BOHR2_TO_M2 = BOHR_TO_M * BOHR_TO_M;
inline constexpr double BOHR2_TO_ANG2 = BOHR_TO_ANG * BOHR_TO_ANG;

// ----------------------------
// Energy conversions
// ----------------------------
inline constexpr double KCAL_PER_KJ = 0.239005736;
inline constexpr double CAL_PER_KJ  = 239.005736;

inline constexpr double KCAL_PER_AU = 627.509474;
inline constexpr double AU_PER_KCAL = 1.0 / KCAL_PER_AU;

// ----------------------------
// Pressure
// ----------------------------
inline constexpr double ATM_TO_PA = 101325.0;

} // namespace pwdft::units

#endif
