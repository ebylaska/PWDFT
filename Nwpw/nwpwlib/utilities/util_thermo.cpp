#include "util_thermo.hpp"
#include <cmath>
#include <iomanip>
#include <array>
#include <algorithm>

/************************************************
 *                                              *
 *          util_molecular_thermochemistry      *
 *                                              *
 ************************************************/
/**
 * @brief Compute molecular thermochemistry from vibrational frequencies.
 *
 * This routine evaluates standard thermochemical quantities for a molecule
 * using the harmonic oscillator approximation applied to a set of vibrational
 * frequencies. The input frequencies are assumed to be obtained from a finite
 * difference (FD) Hessian analysis with translational and rotational modes
 * removed via the Eckart projection.
 *
 * The function computes:
 *  - Zero-point energy (ZPE)
 *  - Thermal correction to internal energy
 *  - Thermal correction to enthalpy
 *  - Translational, rotational, and vibrational entropy
 *  - Constant-volume heat capacity (Cv)
 *
 * Vibrational contributions are evaluated using the standard quantum
 * harmonic oscillator partition function. Translational and rotational
 * contributions are added using classical ideal-gas expressions.
 *
 * The rotational symmetry number sigma is used in the rotational entropy
 * to avoid overcounting indistinguishable molecular orientations.
 *
 * The resulting thermochemical quantities are reported in:
 *  - kcal/mol and Hartree (au) for energies
 *  - cal/mol-K for entropy and heat capacity
 *
 * Frequencies smaller than ~1 cm^-1 are ignored to avoid numerical issues
 * associated with residual rotational or translational modes.
 *
 * @param freq_cm
 *        Vector of vibrational frequencies in cm^-1. These should exclude
 *        the projected translational and rotational modes.
 *
 * @param temperature
 *        Temperature in Kelvin.
 *
 * @param total_mass_au
 *        Molecular mass in atomic units
 *
 * @param pressure
 *        Pressure in Pa.
 *
 * @param sigma
 *        Rotational symmetry number.
 *
 * @param rotor_type
 *        Rotor classification:
 *          0 = atom
 *          1 = linear molecule
 *          2 = nonlinear molecule
 *
 * @param inertia
 *        Principal moments of inertia in amu*Angstrom^2.
 *        For linear molecules, one moment is near zero and the two largest
 *        principal moments are approximately equal.
 *
 * @param out
 *        Output stream used to print formatted thermochemistry results.
 *
 * @ingroup thermochemistry
 */
void util_molecular_thermochemistry(const std::vector<double>& freq_cm,
                                    double temperature,
                                    double total_mass_au,
                                    double pressure,
                                    int sigma,
                                    int rotor_type,
                                    const double inertia[3],
                                    std::ostream& out)
{
    // Thermodynamic constants
    const double R = 8.31446261815324e-3;   // kJ/mol-K
    const double hc_over_k = 1.438776877;   // cm*K

    // SI constants
    const double kB = 1.380649e-23;         // J/K
    const double h  = 6.62607015e-34;       // J*s
    const double pi = 3.14159265358979323846;

    // Unit conversions
    const double me_to_amu   = 1.0 / 1822.888486209;
    const double amu_to_kg   = 1.66053906660e-27;
    const double ang2_to_m2  = 1.0e-20;
    const double kcal_per_kJ = 0.239005736;
    const double cal_per_kJ  = 239.005736;
    const double au_per_kcal = 1.0 / 627.509474;

    // Safety guards
    if (sigma < 1) sigma = 1;
    if (pressure <= 0.0) pressure = 101325.0;
    if (rotor_type < 0 || rotor_type > 2) rotor_type = 2;

    // Copy/sort principal moments for stable handling
    std::array<double,3> I = { inertia[0], inertia[1], inertia[2] };
    std::sort(I.begin(), I.end());

    const bool is_atom   = (rotor_type == 0);
    const bool is_linear = (rotor_type == 1);

    double Ezpe   = 0.0;
    double Evib   = 0.0;   // NWChem-style: includes ZPE
    double Svib   = 0.0;
    double Cv_vib = 0.0;

    // ----------------------------
    // Vibrational contributions
    // ----------------------------
    for (double freq : freq_cm)
    {
        if (freq < 1.0) continue;

        double theta = freq * hc_over_k;   // vibrational temperature (K)
        double x = theta / temperature;

        if (x > 50.0)
        {
            Ezpe += 0.5 * theta;
            Evib += 0.5 * theta;
            continue;
        }

        double ex = std::exp(x);
        double exp_neg_x = std::exp(-x);

        Ezpe += 0.5 * theta;
        Evib += theta * (0.5 + 1.0 / (ex - 1.0));

        Svib += x / (ex - 1.0) - std::log1p(-exp_neg_x);

        Cv_vib += ex * std::pow(x / (ex - 1.0), 2);
    }

    Ezpe   *= R;
    Evib   *= R;
    Svib   *= R;
    Cv_vib *= R;

    // ----------------------------
    // Translational contributions
    // ----------------------------
    double mass_amu = total_mass_au * me_to_amu;
    double mass_kg  = mass_amu * amu_to_kg;
    //double mass_kg = mol_mass * amu_to_kg;


    double qtrans_arg = std::pow(2.0 * pi * mass_kg * kB * temperature / (h * h), 1.5) * (kB * temperature / pressure);

    double Strans   = R * (std::log(qtrans_arg) + 2.5);
    double Etrans   = 1.5 * R * temperature;
    double Cv_trans = 1.5 * R;

    // ----------------------------
    // Rotational contributions
    // ----------------------------
    double Srot   = 0.0;
    double Erot   = 0.0;
    double Cv_rot = 0.0;

    if (is_atom)
    {
        Srot   = 0.0;
        Erot   = 0.0;
        Cv_rot = 0.0;
    }
    else if (is_linear)
    {
        // largest moment after sorting
        double Ilin = I[2];

        if (Ilin > 1.0e-16)
        {
            double I_SI = Ilin * amu_to_kg * ang2_to_m2;
            double theta_rot = h * h / (8.0 * pi * pi * I_SI * kB);

            Srot   = R * (std::log(temperature / (sigma * theta_rot)) + 1.0);
            Erot   = R * temperature;
            Cv_rot = R;
        }
    }
    else
    {

        // Nonlinear molecule
        double IA = I[0] * amu_to_kg * ang2_to_m2;
        double IB = I[1] * amu_to_kg * ang2_to_m2;
        double IC = I[2] * amu_to_kg * ang2_to_m2;

        //if (IA > 1.0e-40 && IB > 1.0e-40 && IC > 1.0e-40)
        if (IA > 0.0 && IB > 0.0 && IC > 0.0)
        {
            double thetaA = h * h / (8.0 * pi * pi * IA * kB);
            double thetaB = h * h / (8.0 * pi * pi * IB * kB);
            double thetaC = h * h / (8.0 * pi * pi * IC * kB);

            double qrot_arg = std::sqrt(pi) * std::pow(temperature, 1.5) / (sigma * std::sqrt(thetaA * thetaB * thetaC));

            Srot   = R * (std::log(qrot_arg) + 1.5);
            Erot   = 1.5 * R * temperature;
            Cv_rot = 1.5 * R;
        }
    }

    // ----------------------------
    // Totals
    // ----------------------------
    double Ethermal = Evib + Etrans + Erot;
    double Hthermal = Ethermal + R * temperature;

    double Stotal   = Svib + Strans + Srot;
    double Cv_total = Cv_vib + Cv_trans + Cv_rot;

    // ----------------------------
    // Output formatting
    // ----------------------------

    out << "\n"
        << " -----------------------------------------\n"
        << "   Canonical Thermochemistry (Ideal Gas)   \n"
        << " -----------------------------------------\n";

    double Ezpe_kcal     = Ezpe     * kcal_per_kJ;
    double Ethermal_kcal = Ethermal * kcal_per_kJ;
    double Hthermal_kcal = Hthermal * kcal_per_kJ;

    double Ezpe_au     = Ezpe_kcal     * au_per_kcal;
    double Ethermal_au = Ethermal_kcal * au_per_kcal;
    double Hthermal_au = Hthermal_kcal * au_per_kcal;

    double Strans_cal = Strans * cal_per_kJ;
    double Srot_cal   = Srot   * cal_per_kJ;
    double Svib_cal   = Svib   * cal_per_kJ;
    double Stotal_cal = Stotal * cal_per_kJ;

    double Cv_trans_cal = Cv_trans * cal_per_kJ;
    double Cv_rot_cal   = Cv_rot   * cal_per_kJ;
    double Cv_vib_cal   = Cv_vib   * cal_per_kJ;
    double Cv_total_cal = Cv_total * cal_per_kJ;

    out << "\n";
    out << " Temperature                        = "
        << std::fixed << std::setprecision(2)
        << temperature << " K\n";

    out << " Pressure                           = "
        << std::fixed << std::setprecision(2)
        << pressure << " Pa\n";

    out << " Rotational symmetry number         = "
        << std::fixed << std::setprecision(0)
        << sigma << "\n";

    out << " Molecule type                      = ";
    if (is_atom)
        out << "atom\n";
    else if (is_linear)
        out << "linear\n";
    else
        out << "nonlinear\n";

    out << std::fixed << std::setprecision(4);
    out << " frequency scaling parameter        =   1.0000\n\n";

    out << std::setprecision(3);
    out << " Zero-Point correction to Energy    = "
        << std::setw(10) << Ezpe_kcal
        << " kcal/mol  ("
        << std::setw(10) << std::setprecision(6) << Ezpe_au
        << " au)\n";

    out << std::setprecision(3);
    out << " Thermal correction to Energy       = "
        << std::setw(10) << Ethermal_kcal
        << " kcal/mol  ("
        << std::setw(10) << std::setprecision(6) << Ethermal_au
        << " au)\n";

    out << std::setprecision(3);
    out << " Thermal correction to Enthalpy     = "
        << std::setw(10) << Hthermal_kcal
        << " kcal/mol  ("
        << std::setw(10) << std::setprecision(6) << Hthermal_au
        << " au)\n";

    out << "\n";

    out << std::setprecision(3);
    out << " Total Entropy                      = "
        << std::setw(10) << Stotal_cal
        << " cal/mol-K\n";

    out << "   - Translational                  = "
        << std::setw(10) << Strans_cal
        << " cal/mol-K\n";

    out << "   - Rotational                     = "
        << std::setw(10) << Srot_cal
        << " cal/mol-K\n";

    out << "   - Vibrational                    = "
        << std::setw(10) << Svib_cal
        << " cal/mol-K\n";

    out << "\n";

    out << " Cv (constant volume heat capacity) = "
        << std::setw(10) << Cv_total_cal
        << " cal/mol-K\n";

    out << "   - Translational                  = "
        << std::setw(10) << Cv_trans_cal
        << " cal/mol-K\n";

    out << "   - Rotational                     = "
        << std::setw(10) << Cv_rot_cal
        << " cal/mol-K\n";

    out << "   - Vibrational                    = "
        << std::setw(10) << Cv_vib_cal
        << " cal/mol-K\n";
}
