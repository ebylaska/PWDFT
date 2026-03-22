#ifndef _UNITS_HPP_
#define _UNITS_HPP_
#pragma once

#pragma once

namespace pwdft::units {

// ----------------------------
// Fundamental constants (SI)
// ----------------------------
static constexpr double KB = 1.380649e-23;          // J/K
static constexpr double H  = 6.62607015e-34;        // J·s
static constexpr double PI = 3.14159265358979323846;
static constexpr double C    = 2.99792458e10;       // cm/s
static constexpr double C_SI = 2.99792458e8;        // m/s
static constexpr double C_CM = 2.99792458e10;       // cm/s

// ----------------------------
// Thermodynamic constants
// ----------------------------
static constexpr double R = 8.31446261815324e-3;    // kJ/mol·K

// ----------------------------
// Spectroscopic conversions
// ----------------------------
static constexpr double HC_OVER_K = (H * C) / KB;   // cm·K
static constexpr double AU_FREQ_TO_CM = 219474.6313705;

// ----------------------------
// Mass conversions
// ----------------------------
static constexpr double AMU_TO_ME = 1822.888486209;
static constexpr double ME_TO_AMU = 1.0 / AMU_TO_ME;
static constexpr double AMU_TO_KG = 1.66053906660e-27;
static constexpr double ME_TO_KG = 9.1093837015e-31;


// ----------------------------
// Length conversions
// ----------------------------
static constexpr double ANG2_TO_M2 = 1.0e-20;
static constexpr double BOHR_TO_M = 5.29177210903e-11;
static constexpr double BOHR_TO_ANG = 0.529177210903;
static constexpr double BOHR2_TO_M2 = BOHR_TO_M * BOHR_TO_M;
static constexpr double BOHR2_TO_ANG2 = BOHR_TO_ANG * BOHR_TO_ANG;

// ----------------------------
// Energy conversions
// ----------------------------
static constexpr double KCAL_PER_KJ = 0.239005736;
static constexpr double CAL_PER_KJ  = 239.005736;

static constexpr double KCAL_PER_AU = 627.509474;
static constexpr double AU_PER_KCAL = 1.0 / KCAL_PER_AU;

// ----------------------------
// Pressure
// ----------------------------
static constexpr double ATM_TO_PA = 101325.0;

}

#endif
