#include "util_thermo.hpp"
#include <cmath>
#include <iomanip>
#include <array>
#include <algorithm>
#include "units.hpp"

using namespace pwdft::units;

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
 *          2 = spherical top
 *          3 = near-spherical top
 *          4 = oblate symmetric top
 *          5 = prolate symmetric top
 *          6 = asymmetric top
 *
 *        Notes:
 *          - Types ≥ 2 are all treated as nonlinear rotors in the
 *            thermodynamic expressions.
 *          - The finer classification (2–6) is used for reporting
 *            and potential symmetry-based validation.
 *
 * @param inertia
 *        Principal moments of inertia in amu*Angstrom^2.
 *        For linear molecules, one moment is near zero and the two largest
 *        principal moments are approximately equal.
 *
 *        Principal moments of inertia (I_A, I_B, I_C) in amu·Å²,
 *        corresponding to the principal axes obtained from the
 *        inertia tensor diagonalization.
 *
 *        The ordering is preserved as provided (no sorting is applied),
 *        so the associated rotational constants (A, B, C) reflect the
 *        geometry-defined principal axes.
 *
 *        Interpretation:
 *        - Atom: all moments ≈ 0
 *        - Linear: one moment ≈ 0, two large moments
 *        - Nonlinear: all moments > 0
 *
 * @param out
 *        Output stream used to print formatted thermochemistry results.
 *
 * @ingroup thermochemistry
 */
ThermoResults util_molecular_thermochemistry(const std::vector<double>& freq_cm,
                                             double temperature,
                                             double total_mass_au,
                                             double pressure,
                                             int sigma,
                                             int rotor_type,
                                             const double inertia[3],
                                             std::ostream* out)
{

    // Safety guards
    if (sigma < 1) sigma = 1;
    //if (pressure <= 0.0) pressure = 101325.0;
    if (pressure <= 0.0) pressure = ATM_TO_PA;
    //if (rotor_type < 0 || rotor_type > 2) rotor_type = 2;

    // Copy/sort principal moments for stable handling
    std::array<double,3> I = { inertia[0], inertia[1], inertia[2] };
    //std::sort(I.begin(), I.end());

    const bool is_atom   = (rotor_type == 0);
    const bool is_linear = (rotor_type == 1);

    std::string rotor_class;
    if (rotor_type==0) rotor_class = "atom";
    if (rotor_type==1) rotor_class = "linear rotor";
    if (rotor_type==2) rotor_class = "spherical top";
    if (rotor_type==3) rotor_class = "near-spherical top";
    if (rotor_type==4) rotor_class = "oblate symmetric top";
    if (rotor_type==5) rotor_class = "prolate symmetric top";
    if (rotor_type==6) rotor_class = "asymmetric top";


    double Ezpe   = 0.0;
    double Evib   = 0.0;   // NWChem-style: includes ZPE
    double Svib   = 0.0;
    double Cv_vib = 0.0;

    // ----------------------------
    // Vibrational contributions
    // ----------------------------
    for (double freq : freq_cm)
    {
        //if (freq < 1.0) continue;
        if (freq < 1e-3) continue;

        //double theta = freq * hc_over_k;   // vibrational temperature (K)
        double theta = freq * HC_OVER_K;
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
    double mass_amu = total_mass_au * ME_TO_AMU;
    double mass_kg  = mass_amu * AMU_TO_KG;

    //double qtrans_arg = std::pow(2.0 * pi * mass_kg * kB * temperature / (h * h), 1.5) * (kB * temperature / pressure);
    double qtrans_arg = std::pow(2.0 * PI * mass_kg * KB * temperature / (H * H), 1.5) * (KB * temperature / pressure);

    double Strans   = R * (std::log(qtrans_arg) + 2.5);
    double Etrans   = 1.5 * R * temperature;
    double Cv_trans = 1.5 * R;

    // ----------------------------------------------
    // Rotational contributions (stable formulation)
    // ----------------------------------------------
    double Srot   = 0.0;
    double Erot   = 0.0;
    double Cv_rot = 0.0;

    const double C_inertia = AMU_TO_KG * ANG2_TO_M2;

    const double I_tol = 1e-12; // in amu·Å² space
    const double rot_cm_prefactor = H / (8.0 * PI * PI * C_CM * C_inertia);   // C must be cm/s
    const double rot_K_prefactor  = H * H / (8.0 * PI * PI * KB * C_inertia);

    //double IA = I[0] * C_inertia;
    //double IB = I[1] * C_inertia;
    //double IC = I[2] * C_inertia;

    double A_cm = -99.9;
    double B_cm = -99.9;
    double C_cm = -99.9;
    double A_K  = -99.9;
    double B_K  = -99.9;
    double C_K  = -99.9;
    if (I[0] > I_tol) 
    {
       A_cm = rot_cm_prefactor / I[0];
       A_K  = rot_K_prefactor  / I[0];
    }
    if (I[1] > I_tol)  
    {
       B_cm = rot_cm_prefactor / I[1];
       B_K  = rot_K_prefactor  / I[1];
    }
    if (I[2] > I_tol)  
    {
       C_cm = rot_cm_prefactor / I[2];
       C_K  = rot_K_prefactor  / I[2];
    }



    if (is_atom)
    {
       Srot   = 0.0;
       Erot   = 0.0;
       Cv_rot = 0.0;
    }
    else if (is_linear)
    {
       // largest principal moment for a linear rotor
       //double Ilin = I[2];
       double Ilin = std::max({I[0], I[1], I[2]});

       if (Ilin > I_tol)
       {
           const double log_rot_prefactor_linear = std::log(8.0 * PI * PI * KB / (H * H)) + std::log(C_inertia);
           //double I_SI = Ilin * C_inertia;
 
           // qrot = T / (sigma * theta_rot)
           //      = (8*pi^2*kB*T*I)/(sigma*h^2)
           //double log_qrot =
           //      std::log(temperature)
           //    - std::log((double)sigma)
           //    + std::log(8.0 * PI * PI * KB / (H * H))
           //    + std::log(I_SI);
 
           double log_qrot = log_rot_prefactor_linear
                           - std::log((double)sigma)
                           + std::log(temperature)
                           + std::log(Ilin);
 
           Srot   = R * (log_qrot + 1.0);
           Erot   = R * temperature;
           Cv_rot = R;
       }
    }

    // Nonlinear molecule
    else
    {
       if (I[0] > I_tol && I[1] > I_tol && I[2] > I_tol)
       {
           // qrot = sqrt(pi)/sigma * (8*pi^2*kB*T/h^2)^(3/2) * sqrt(IA*IB*IC)
           //double log_qrot =
           //      0.5 * std::log(PI)
           //    - std::log((double)sigma)
           //    + 1.5 * std::log(8.0 * PI * PI * KB / (H * H))
           //    + 1.5 * std::log(temperature)
           //    + 0.5 * std::log(IA * IB * IC);
 
           double log_qrot = 0.5 * std::log(PI)
                           - std::log((double)sigma)
                           + 1.5 * std::log(8.0 * PI * PI * KB / (H * H))
                           + 1.5 * std::log(temperature)
                           + 0.5 * (std::log(I[0]) + std::log(I[1]) + std::log(I[2]))
                           + 1.5 * std::log(C_inertia);
 
           Srot   = R * (log_qrot + 1.5);
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

    double Gthermal = Hthermal - temperature * Stotal;

    double Ezpe_kcal     = Ezpe     * KCAL_PER_KJ;
    double Ethermal_kcal = Ethermal * KCAL_PER_KJ;
    double Hthermal_kcal = Hthermal * KCAL_PER_KJ;
    double Gthermal_kcal = Gthermal * KCAL_PER_KJ;

    double Ezpe_au     = Ezpe_kcal     * AU_PER_KCAL;
    double Ethermal_au = Ethermal_kcal * AU_PER_KCAL;
    double Hthermal_au = Hthermal_kcal * AU_PER_KCAL;
    double Gthermal_au = Gthermal_kcal * AU_PER_KCAL;

    double Strans_cal = Strans * CAL_PER_KJ;
    double Srot_cal   = Srot   * CAL_PER_KJ;
    double Svib_cal   = Svib   * CAL_PER_KJ;
    double Stotal_cal = Stotal * CAL_PER_KJ;

    double Cv_trans_cal = Cv_trans * CAL_PER_KJ;
    double Cv_rot_cal   = Cv_rot   * CAL_PER_KJ;
    double Cv_vib_cal   = Cv_vib   * CAL_PER_KJ;
    double Cv_total_cal = Cv_total * CAL_PER_KJ;

    // ----------------------------
    // Output formatting
    // ----------------------------
    if (out)
    {
       auto& myout = *out;

        myout << "\n"
              << " -----------------------------------------\n"
              << "   Canonical Thermochemistry (Ideal Gas)   \n"
              << " -----------------------------------------\n";
 
 
        myout << "\n";
        myout << " temperature                        = "
              << std::fixed << std::setprecision(2)
              << temperature << " K\n";
 
        myout << " pressure                           = "
              << std::fixed << std::setprecision(2)
              << pressure << " Pa (" 
              << pressure/ATM_TO_PA << " atm)\n";
 
        myout << " molecular mass                     = "
              << std::fixed << std::setprecision(3)
              << mass_amu << " amu\n";
 
        myout << " molecular mass (au, me)            = "
              << std::fixed << std::setprecision(3)
              << total_mass_au << "\n";
 
        myout << " rotational symmetry number         = "
              << std::fixed << std::setprecision(0)
              << sigma << "\n";
 
        myout << " rotor type                         = ";
        if (is_atom)
        {
            myout << "atom\n";
            myout << " (atom: no rotational constants)\n";
        }
        else if (is_linear)
            myout << "linear\n";
        else
            myout << "nonlinear\n";
 
        myout << " rotor class                        = " << rotor_class << "\n";
        myout << " Rotational Constants\n";
        myout << std::fixed << std::setprecision(6);
 
        if ((A_cm>0.0) && (A_K>0.0))
           myout << "   - A = " << std::setw(10) << A_cm << " cm-1  (" << std::setw(10) << A_K << " K)\n";
        if ((B_cm>0.0) && (B_K>0.0))
           myout << "   - B = " << std::setw(10) << B_cm << " cm-1  (" << std::setw(10) << B_K << " K)\n";
        if ((C_cm>0.0) && (C_K>0.0))
           myout << "   - C = " << std::setw(10) << C_cm << " cm-1  (" << std::setw(10) << C_K << " K)\n";
 
        myout << std::fixed << std::setprecision(4);
        myout << " frequency scaling parameter        =   1.0000\n\n";
 
        myout << std::setprecision(3);
        myout << " zero-point correction to energy         = "
              << std::setw(10) << Ezpe_kcal
              << " kcal/mol  ("
              << std::setw(10) << std::setprecision(6) << Ezpe_au
              << " au)\n";
 
        myout << std::setprecision(3);
        myout << " thermal correction to energy            = "
              << std::setw(10) << Ethermal_kcal
              << " kcal/mol  ("
              << std::setw(10) << std::setprecision(6) << Ethermal_au
              << " au)\n";
 
        myout << std::setprecision(3);
        myout << " thermal correction to enthalpy          = "
              << std::setw(10) << Hthermal_kcal
              << " kcal/mol  ("
              << std::setw(10) << std::setprecision(6) << Hthermal_au
              << " au)\n";
 
        myout << std::setprecision(3);
        myout << " thermal correction to Gibbs free energy = "
              << std::setw(10) << Gthermal_kcal
              << " kcal/mol  ("
              << std::setw(10) << std::setprecision(6) << Gthermal_au
              << " au)\n";
 
 
        myout << "\n";
 
        myout << std::setprecision(3);
        myout << " total entropy                           = "
              << std::setw(10) << Stotal_cal
              << " cal/mol-K\n";
 
        myout << "   - translational                       = "
              << std::setw(10) << Strans_cal
              << " cal/mol-K\n";
 
        myout << "   - rotational                          = "
              << std::setw(10) << Srot_cal
              << " cal/mol-K\n";
 
        myout << "   - vibrational                         = "
              << std::setw(10) << Svib_cal
              << " cal/mol-K\n";
 
        myout << "\n";
 
        myout << " Cv (constant volume heat capacity)      = "
              << std::setw(10) << Cv_total_cal
              << " cal/mol-K\n";
 
        myout << "   - translational                       = "
              << std::setw(10) << Cv_trans_cal
              << " cal/mol-K\n";
 
        myout << "   - rotational                          = "
              << std::setw(10) << Cv_rot_cal
              << " cal/mol-K\n";
 
        myout << "   - vibrational                         = "
              << std::setw(10) << Cv_vib_cal
              << " cal/mol-K\n";
    }

    ThermoResults result{};

    result.Ezpe_au      = Ezpe_au;
    result.Ethermal_au  = Ethermal_au;
    result.Hthermal_au  = Hthermal_au;
    result.Gthermal_au  = Gthermal_au;

    result.Ezpe_kcal     = Ezpe_kcal;
    result.Ethermal_kcal = Ethermal_kcal;
    result.Hthermal_kcal = Hthermal_kcal;
    result.Gthermal_kcal = Gthermal_kcal;

    result.S_total = Stotal_cal;
    result.S_trans = Strans_cal;
    result.S_rot   = Srot_cal;
    result.S_vib   = Svib_cal;

    result.Cv_total = Cv_total_cal;
    result.Cv_trans = Cv_trans_cal;
    result.Cv_rot   = Cv_rot_cal;
    result.Cv_vib   = Cv_vib_cal;

    result.A_cm = A_cm;
    result.B_cm = B_cm;
    result.C_cm = C_cm;

    result.sigma      = sigma;
    result.rotor_type = rotor_type;

    result.temperature = temperature;
    result.pressure    = pressure;

    return result;
}
