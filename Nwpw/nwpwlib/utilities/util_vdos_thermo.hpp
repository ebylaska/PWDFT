#ifndef UTIL_VDOS_THERMO_HPP
#define UTIL_VDOS_THERMO_HPP

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>
#include "util_ascii_plot.hpp"
#include "util_thermo.hpp"


namespace pwdft {

/******************************************
 *                                        *
 *            _util_vdos_thermo_          *
 *                                        *
 ******************************************/
/*
 * Compute thermodynamic functions from the vibrational density of states.
 *
 * Input:
 *   efn is a flat array storing npoints x 3 data:
 *
 *     efn[3*k + 0] = w_k    frequency in atomic units
 *     efn[3*k + 1] = g(w_k) vibrational density of states
 *     efn[3*k + 2] = N(w_k) integrated number of vibrational states
 *
 * Parameters:
 *   scalefreq      - vibrational scaling factor
 *   threshfreq_cm  - threshold frequency in cm^-1 for splitting low/high
 *   temperatures   - temperatures in K
 *
 * Notes:
 *   - efn is temporarily rescaled for plotting and then restored.
 *   - Frequencies are assumed to be in atomic units on input.
 */
inline std::vector<ThermoResults> util_vdos_thermo(int nion,
                                                   int npoints,
                                                   double* efn,
                                                   double scalefreq,
                                                   double threshfreq_cm,
                                                   const std::vector<double>& temperatures,
                                                   std::ostream& coutput)
{

    (void)nion;

    constexpr double AUKCAL = 627.5093314;
    constexpr double Rgas   = 1.9863e0 / 1000.0 / AUKCAL;
    constexpr double autocm = 219474.6313705;

    const double threshfreq = threshfreq_cm / autocm;

    auto E = [efn](int k) -> double& { return efn[3 * k + 0]; };
    auto G = [efn](int k) -> double& { return efn[3 * k + 1]; };
    auto N = [efn](int k) -> double& { return efn[3 * k + 2]; };

    std::vector<ThermoResults> results;
    results.reserve(temperatures.size());

    /************************************************************
     *************** ascii plot of the density of states ********
     ************************************************************/
    coutput << "\n\n\n";
    coutput << "     g(w) and N(w) Used For Calculating Thermodynamic Functions\n";
    coutput << "     ----------------------------------------------------------\n";

    double emin =  9.99999999e6;
    double gmin =  9.99999999e6;
    double nmin =  9.99999999e6;
    double emax = -9.99999999e6;
    double gmax = -9.99999999e6;
    double nmax = -9.99999999e6;

    for (int k = 0; k < npoints; ++k) {
        const double xx = E(k);
        const double yy = G(k);
        const double zz = N(k);

        emin = std::min(emin, xx);
        gmin = std::min(gmin, yy);
        nmin = std::min(nmin, zz);

        emax = std::max(emax, xx);
        gmax = std::max(gmax, yy);
        nmax = std::max(nmax, zz);
    }

    coutput << "\n      Density of States Details\n";
    coutput << "         - grid size (npoints)         = "
            << std::setw(10) << npoints << "\n";
    coutput << "         - minimum frequency w         = "
            << std::fixed << std::setprecision(1)
            << std::setw(10) << emin * autocm << " cm-1\n";
    coutput << "         - maximum frequency w         = "
            << std::fixed << std::setprecision(1)
            << std::setw(10) << emax * autocm << " cm-1\n";
    coutput << "         - minimum g(w)                = "
            << std::scientific << std::setprecision(3)
            << std::setw(10) << gmin / autocm << " states/cm-1\n";
    coutput << "         - maximum g(w)                = "
            << std::scientific << std::setprecision(3)
            << std::setw(10) << gmax / autocm << " states/cm-1\n";
    coutput << "         - minimum N(w)                = "
            << std::fixed << std::setprecision(1)
            << std::setw(10) << nmin << " states\n";
    coutput << "         - maximum N(w)                = "
            << std::fixed << std::setprecision(1)
            << std::setw(10) << nmax << " states\n\n\n";

    for (int k = 0; k < npoints; ++k) {
        E(k) *= autocm;
        G(k) /= autocm;
    }

    if (std::abs(gmax - gmin) < 1.0e-3)
        gmax = gmin + 1.0e-3;

    {
        std::vector<double> x(npoints), y(npoints);
        for (int k = 0; k < npoints; ++k) {
            x[k] = E(k);
            y[k] = G(k);
        }

        util_ascii::util_ascii_setwindow(emin * autocm, emax * autocm,
                                         gmin / autocm, gmax / autocm);

        util_ascii::util_ascii_plotter(" ", coutput, x, y, "*",
                                       "g(w) - Vibrational Density of States",
                                       "w (cm-1)",
                                       "g(w) (states/cm-1)");
    }

    coutput << "\n\n";

    {
        std::vector<double> x(npoints), y(npoints);
        for (int k = 0; k < npoints; ++k) {
            x[k] = E(k);
            y[k] = N(k);
        }

        util_ascii::util_ascii_setwindow(emin * autocm, emax * autocm,
                                         nmin, nmax);

        util_ascii::util_ascii_plotter(" ", coutput, x, y, "*",
                                       "N(w) - Number of Vibrational States",
                                       "w (cm-1)",
                                       "N(w) (states)");
    }

    coutput << "\n";

    for (int k = 0; k < npoints; ++k) {
        E(k) /= autocm;
        G(k) *= autocm;
    }

    /************************************************************
     *************** Compute the Thermodynamic Functions ********
     ************************************************************/
    coutput << "\n\n";
    coutput << "     Thermodynamic Functions from Vibrational Density of States\n";
    coutput << "     ----------------------------------------------------------\n";

    for (double temp : temperatures) {
        const double RT  = Rgas * temp;
        const double RT2 = 1.0 / (2.0 * RT);
        const double RT1 = 1.0 / RT;

        double ezero       = 0.0;
        double ezero1      = 0.0;
        double ezero2      = 0.0;
        double ethermal    = 0.0;
        double ethermal1   = 0.0;
        double ethermal2   = 0.0;
        double fhelmholtz  = 0.0;
        double fhelmholtz1 = 0.0;
        double fhelmholtz2 = 0.0;
        double sa          = 0.0;
        double sa1         = 0.0;
        double sa2         = 0.0;
        double sb          = 0.0;
        double sb1         = 0.0;
        double sb2         = 0.0;
        double Svib        = 0.0;
        double Svib1       = 0.0;
        double Svib2       = 0.0;
        double Cv_vib      = 0.0;
        double Cv_vib1     = 0.0;
        double Cv_vib2     = 0.0;

        const double de = E(1) - E(0);

        for (int k = 1; k < npoints - 1; ++k) {
            const double hbarw = scalefreq * E(k);
            const double g = 0.5 * (N(k + 1) - N(k - 1)) / de;

            double xdum = (temp > 0.0) ? std::exp(-hbarw * RT1) : 0.0;
            xdum = xdum / (1.0 - xdum);

            ezero      += g * hbarw * 0.5;
            ethermal   += g * hbarw * (0.5 + xdum);
            fhelmholtz += g * RT * std::log(2.0 * std::sinh(hbarw * RT2));

            sa     += g * hbarw * xdum;
            sb     += g * std::log(1.0 - std::exp(-hbarw * RT1));
            Cv_vib += g * std::exp(hbarw * RT1) * std::pow(hbarw * RT1 * xdum, 2);

            if (hbarw < threshfreq) {
                ezero1      += g * hbarw * 0.5;
                ethermal1   += g * hbarw * (0.5 + xdum);
                fhelmholtz1 += g * RT * std::log(2.0 * std::sinh(hbarw * RT2));
                sa1         += g * hbarw * xdum;
                sb1         += g * std::log(1.0 - std::exp(-hbarw * RT1));
                Cv_vib1     += g * std::exp(hbarw * RT1) * std::pow(hbarw * RT1 * xdum, 2);
            } else {
                ezero2      += g * hbarw * 0.5;
                ethermal2   += g * hbarw * (0.5 + xdum);
                fhelmholtz2 += g * RT * std::log(2.0 * std::sinh(hbarw * RT2));
                sa2         += g * hbarw * xdum;
                sb2         += g * std::log(1.0 - std::exp(-hbarw * RT1));
                Cv_vib2     += g * std::exp(hbarw * RT1) * std::pow(hbarw * RT1 * xdum, 2);
            }
        }

        ezero       *= de;
        ezero1      *= de;
        ezero2      *= de;
        ethermal    *= de;
        ethermal1   *= de;
        ethermal2   *= de;
        fhelmholtz  *= de;
        fhelmholtz1 *= de;
        fhelmholtz2 *= de;

        if (temp > 0.0) {
            Svib  = sa  * de / temp - Rgas * sb  * de;
            Svib1 = sa1 * de / temp - Rgas * sb1 * de;
            Svib2 = sa2 * de / temp - Rgas * sb2 * de;
        }

        Cv_vib  = Rgas * Cv_vib  * de;
        Cv_vib1 = Rgas * Cv_vib1 * de;
        Cv_vib2 = Rgas * Cv_vib2 * de;

        ThermoResults res{};

        res.Ezpe_au      = ezero;
        res.Ethermal_au  = ethermal;
        res.Hthermal_au  = ethermal;
        res.Gthermal_au  = fhelmholtz;

        res.Ezpe_kcal      = ezero * AUKCAL;
        res.Ethermal_kcal  = ethermal * AUKCAL;
        res.Hthermal_kcal  = ethermal * AUKCAL;
        res.Gthermal_kcal  = fhelmholtz * AUKCAL;

        res.S_vib   = Svib * AUKCAL * 1000.0;
        res.S_total = res.S_vib;
        res.S_trans = 0.0;
        res.S_rot   = 0.0;

        res.Cv_vib   = Cv_vib * AUKCAL * 1000.0;
        res.Cv_total = res.Cv_vib;
        res.Cv_trans = 0.0;
        res.Cv_rot   = 0.0;

        res.A_cm = 0.0;
        res.B_cm = 0.0;
        res.C_cm = 0.0;

        res.sigma = 0;
        res.rotor_type = 0;

        res.temperature = temp;
        res.pressure = 0.0;
        results.push_back(res);

        coutput << "\n      Temperature                      = "
                << std::fixed << std::setprecision(2) << temp
                << "K (RT=" << std::scientific << std::setprecision(6)
                << RT << " au)\n";

        coutput << "      frequency scaling parameter      = "
                << std::fixed << std::setprecision(4) << scalefreq << "\n";

        coutput << "      frequency dividing parameter     = "
                << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1\n";

        coutput << "\n      Helmholtz Free Energy            = "
                << std::fixed << std::setprecision(3) << fhelmholtz * AUKCAL
                << " kcal/mol  (" << std::setprecision(6) << fhelmholtz << " au)\n";

        coutput << "        - (less than   " << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)  = " << std::setprecision(3)
                << fhelmholtz1 * AUKCAL << " kcal/mol  (" << std::setprecision(6)
                << fhelmholtz1 << " au)\n";

        coutput << "        - (greater than" << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)  = " << std::setprecision(3)
                << fhelmholtz2 * AUKCAL << " kcal/mol  (" << std::setprecision(6)
                << fhelmholtz2 << " au)\n";

        coutput << "\n      Zero-Point correction to Energy  = "
                << std::fixed << std::setprecision(3) << ezero * AUKCAL
                << " kcal/mol  (" << std::setprecision(6) << ezero << " au)\n";

        coutput << "        - (less than   " << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)  = " << std::setprecision(3)
                << ezero1 * AUKCAL << " kcal/mol  (" << std::setprecision(6)
                << ezero1 << " au)\n";

        coutput << "        - (greater than" << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)  = " << std::setprecision(3)
                << ezero2 * AUKCAL << " kcal/mol  (" << std::setprecision(6)
                << ezero2 << " au)\n";

        coutput << "\n      Thermal correction to Energy     = "
                << std::fixed << std::setprecision(3) << ethermal * AUKCAL
                << " kcal/mol  (" << std::setprecision(6) << ethermal << " au)\n";

        coutput << "        - (less than   " << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)  = " << std::setprecision(3)
                << ethermal1 * AUKCAL << " kcal/mol  (" << std::setprecision(6)
                << ethermal1 << " au)\n";

        coutput << "        - (greater than" << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)  = " << std::setprecision(3)
                << ethermal2 * AUKCAL << " kcal/mol  (" << std::setprecision(6)
                << ethermal2 << " au)\n";

        coutput << "\n      Entropy                            = "
                << std::fixed << std::setprecision(3)
                << Svib * AUKCAL * 1000.0 << " cal/mol-K\n";

        coutput << "        - (less than   " << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)    = " << std::setprecision(3)
                << Svib1 * AUKCAL * 1000.0 << " cal/mol-K\n";

        coutput << "        - (greater than" << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)    = " << std::setprecision(3)
                << Svib2 * AUKCAL * 1000.0 << " cal/mol-K\n";

        coutput << "\n      Cv (constant volume heat capacity) = "
                << std::fixed << std::setprecision(3)
                << Cv_vib * AUKCAL * 1000.0 << " cal/mol-K\n";

        coutput << "        - (less than   " << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)    = " << std::setprecision(3)
                << Cv_vib1 * AUKCAL * 1000.0 << " cal/mol-K\n";

        coutput << "        - (greater than" << std::fixed << std::setprecision(1)
                << threshfreq_cm << " cm-1)    = " << std::setprecision(3)
                << Cv_vib2 * AUKCAL * 1000.0 << " cal/mol-K\n\n";
    }

    coutput << std::defaultfloat;
 
    return results;
}

} // namespace pwdft

#endif
