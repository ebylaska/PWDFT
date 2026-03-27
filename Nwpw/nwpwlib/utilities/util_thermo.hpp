#ifndef _UTIL_THERMO_HPP_
#define _UTIL_THERMO_HPP_

#pragma once

#include <vector>
#include <iostream>

struct ThermoResults
{
    // energies (au)
    double Ezpe_au;
    double Ethermal_au;
    double Hthermal_au;
    double Gthermal_au;

    // energies (kcal/mol)
    double Ezpe_kcal;
    double Ethermal_kcal;
    double Hthermal_kcal;
    double Gthermal_kcal;

    // entropy (cal/mol-K)
    double S_total;
    double S_trans;
    double S_rot;
    double S_vib;

    // heat capacity (cal/mol-K)
    double Cv_total;
    double Cv_trans;
    double Cv_rot;
    double Cv_vib;

    // rotational constants (cm^-1)
    double A_cm, B_cm, C_cm;

    // metadata
    int sigma;
    int rotor_type;

    // thermodynamic conditions
    double temperature; // K
    double pressure;    // Pa
};

//void util_molecular_thermochemistry(const std::vector<double>&, double, double, std::ostream&);

//void util_molecular_thermochemistry(const std::vector<double>&, double, double, double, int, int, const double *, std::ostream&);

ThermoResults util_molecular_thermochemistry(const std::vector<double>& freq_cm, 
                                             double temperature, 
                                             double mol_mass, 
                                             double pressure, 
                                             int sigma, 
                                             int rotor_type, 
                                             const double* inertia,
                                             std::ostream* out = nullptr);


#endif
