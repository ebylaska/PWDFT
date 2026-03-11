#include "util_thermo.hpp"
#include <cmath>
#include <iomanip>

/************************************************
 *                                              *
 *          util_molecular_thermochemistry      *
 *                                              *
 ************************************************/

void util_molecular_thermochemistry(const std::vector<double>& freq_cm,
                                    double temperature,
                                    double mol_mass,
                                    std::ostream& out)
{
    const double R = 8.31446261815324e-3; // kJ/mol-K
    const double hc_over_k = 1.438776877; // cm*K

    const double kcal_per_kJ = 0.239005736; // kJ → kcal
    const double cal_per_kJ  = 239.005736;  // kJ → cal
    const double au_per_kcal = 1.0 / 627.509474;

    double Ezpe = 0.0;
    double Evib = 0.0;
    double Svib = 0.0;
    double Cv_vib = 0.0;

    for(double freq : freq_cm)
    {
        if(freq < 1.0) continue;

        double theta = freq * hc_over_k;
        double x = theta/temperature;

        if(x > 50.0)
        {
            Ezpe += 0.5 * theta;
            Evib += 0.5 * theta;
            continue;
        }

        double ex = std::exp(x);

        Ezpe += 0.5 * theta;
        Evib += theta*(0.5 + 1.0/(ex-1.0));

        Svib += x/(ex-1.0) - std::log(1.0-std::exp(-x));

        Cv_vib += ex * std::pow(x/(ex-1.0),2);
    }

    Ezpe *= R;
    Evib *= R;

    Svib *= R;
    Cv_vib *= R;

    double Etrans = 1.5*R*temperature;
    double Erot   = 1.5*R*temperature;

    double Cv_trans = 1.5*R;
    double Cv_rot   = 1.5*R;

    double Ethermal = Evib + Etrans + Erot;
    double Hthermal = Ethermal + R*temperature;

    double Cv_total = Cv_vib + Cv_trans + Cv_rot;

    //out << "\n Thermochemistry (T=" << temperature << " K)\n";
    //out << " -----------------------------------------\n";


    out << "\n Vibrational analysis via the FX method" << std::endl; 
    out << "--- with translations and rotations projected out ---" << std::endl;
    out << "--- via the Eckart algorithm                      ---" << std::endl;        
    out << " -----------------------------------------\n";

    double Ezpe_kcal = Ezpe * kcal_per_kJ;
    double Ethermal_kcal = Ethermal * kcal_per_kJ;
    double Hthermal_kcal = Hthermal * kcal_per_kJ;

    double Ezpe_au = Ezpe_kcal * au_per_kcal;
    double Ethermal_au = Ethermal_kcal * au_per_kcal;
    double Hthermal_au = Hthermal_kcal * au_per_kcal;

    double Svib_cal = Svib * cal_per_kJ;
    double Cv_total_cal = Cv_total * cal_per_kJ;

    out << "\n";
    out << " Temperature                        = "
        << std::fixed << std::setprecision(2)
        << temperature << " K\n";

    out << " frequency scaling parameter        =   1.0000\n\n";

    out << " Zero-Point correction to Energy    = "
        << std::setw(10) << std::setprecision(3) << Ezpe_kcal
        << " kcal/mol  ("
        << std::setw(10) << std::setprecision(6) << Ezpe_au
        << " au)\n";

    out << " Thermal correction to Energy       = "
        << std::setw(10) << std::setprecision(3) << Ethermal_kcal
        << " kcal/mol  ("
        << std::setw(10) << std::setprecision(6) << Ethermal_au
        << " au)\n";

    out << " Thermal correction to Enthalpy     = "
        << std::setw(10) << std::setprecision(3) << Hthermal_kcal
        << " kcal/mol  ("
        << std::setw(10) << std::setprecision(6) << Hthermal_au
        << " au)\n";

    out << "\n";

    out << " Total Entropy                      = "
        << std::setw(10) << std::setprecision(3)
        << Svib_cal
        << " cal/mol-K\n";

    out << "\n";

    out << " Cv (constant volume heat capacity) = "
        << std::setw(10) << std::setprecision(3)
        << Cv_total_cal
        << " cal/mol-K\n";

    /*
    out << " Zero-point energy        = " << Ezpe << " kJ/mol\n";
    out << " Thermal energy correction= " << Ethermal << " kJ/mol\n";
    out << " Enthalpy correction      = " << Hthermal << " kJ/mol\n";

    out << "\n Heat capacity Cv         = " << Cv_total << " kJ/mol-K\n";
    out << " Vibrational entropy Svib = " << Svib << " kJ/mol-K\n";
    */
}
