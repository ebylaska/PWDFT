#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <map>

#include "json.hpp"
#include "parsestring.hpp"

using json = nlohmann::json;
// using ordered_json = nlohmann::ordered_json;
// using json = nlohmann::ordered_json;

namespace pwdft {

json periodic_table_mass = json::parse(
    "{ \"H\"  : 1.008, \"He\" : 4.0026, \"Li\" : 7.016, \"Be\" : 9.01218, "
    "\"B\"  : 11.00931, \"C\"  : 12.0, \"N\"  : 14.00307, \"O\"  : 15.99491, "
    "\"F\"  : 18.9984, \"Ne\" : 19.99244, \"Na\" : 22.9898, \"Mg\" : 23.98504, "
    "\"Al\" : 26.98154, \"Si\" : 27.97693, \"P\"  : 30.97376, \"S\"  : "
    "31.97207, \"Cl\" : 34.96885, \"Ar\" : 39.9624, \"K\"  : 38.96371, \"Ca\" "
    ": 39.96259, \"Sc\" : 44.95592, \"Ti\" : 45.948, \"V\"  : 50.9440, \"Cr\" "
    ": 51.9405, \"Mn\" : 54.9381, \"Fe\" : 55.9349, \"Co\" : 58.9332, \"Ni\" : "
    "57.9353, \"Cu\" : 62.9298, \"Zn\" : 63.9291, \"Ga\" : 68.9257, \"Ge\" : "
    "73.9219, \"As\" : 74.9216, \"Se\" : 78.9183, \"Br\" : 79.9165, \"Kr\" : "
    "83.912, \"Rb\" : 84.9117, \"Sr\" : 87.9056, \"Y\"  : 88.9054, \"Zr\" : "
    "89.9043, \"Nb\" : 92.9060, \"Mo\" : 97.9055, \"Tc\" : 97.9072, \"Ru\" : "
    "101.9037, \"Rh\" : 102.9048, \"Pd\" : 105.9032, \"Ag\" : 106.90509, "
    "\"Cd\" : 113.9036, \"In\" : 114.9041, \"Sn\" : 117.9018, \"Sb\" : "
    "120.9038, \"Te\" : 129.9067, \"I\"  : 126.9004, \"Xe\" : 131.9042, \"Cs\" "
    ": 132.9051, \"Ba\" : 137.9050, \"La\" : 138.9061, \"Ce\" : 139.9053, "
    "\"Pr\" : 140.9074, \"Nd\" : 143.9099, \"Pm\" : 144.9128, \"Sm\" : "
    "151.9195, \"Eu\" : 152.920, \"Gd\" : 157.9241, \"Tb\" : 159.9250, \"Dy\" "
    ": 163.9288, \"Ho\" : 164.9303, \"Er\" : 165.930, \"Tm\" : 168.9344, "
    "\"Yb\" : 173.9390, \"Lu\" : 174.9409, \"Hf\" : 179.9468, \"Ta\" : "
    "180.948, \"W\"  : 183.9510, \"Re\" : 186.9560, \"Os\" : 189.9586, \"Ir\" "
    ": 192.9633, \"Pt\" : 194.9648, \"Au\" : 196.9666, \"Hg\" : 201.9706, "
    "\"Tl\" : 204.9745, \"Pb\" : 207.9766, \"Bi\" : 208.9804, \"Po\" : "
    "209.9829, \"At\" : 210.9875, \"Rn\" : 222.0175, \"Fr\" : 223.0198, \"Ra\" "
    ": 226.0254, \"Ac\" : 227.0278, \"Th\" : 232.0382, \"Pa\" : 231.0359, "
    "\"U\"  : 238.0508, \"Np\" : 237.0482, \"Pu\" : 244.0642, \"Am\" : "
    "243.0614, \"Cm\" : 247.0704, \"Bk\" : 247.0703, \"Cf\" : 251.0796, \"Es\" "
    ": 252.0829, \"Fm\" : 257.0950, \"Md\" : 258.0986, \"No\" : 259.1009, "
    "\"Lr\" : 262.1100, \"Rf\" : 261.1087, \"Ha\" : 262.1138, \"Sg\" : "
    "266.1219, \"Bh\" : 262.1229, \"Hs\" : 267.1318, \"Mt\" : 268.1388 }");

json periodic_table_Z = json::parse(
    "{ \"H\"  : 1, \"He\" : 2, \"Li\" : 3, \"Be\" : 4, \"B\"  : 5, \"C\"  : 6, "
    "\"N\"  : 7, \"O\"  : 8, \"F\"  : 9, \"Ne\" : 10, \"Na\" : 11, \"Mg\" : "
    "12, \"Al\" : 13, \"Si\" : 14, \"P\"  : 15, \"S\"  : 16, \"Cl\" : 17, "
    "\"Ar\" : 18, \"K\"  : 19, \"Ca\" : 20, \"Sc\" : 21, \"Ti\" : 22, \"V\"  : "
    "23, \"Cr\" : 24, \"Mn\" : 25, \"Fe\" : 26, \"Co\" : 27, \"Ni\" : 28, "
    "\"Cu\" : 29, \"Zn\" : 30, \"Ga\" : 31, \"Ge\" : 32, \"As\" : 33, \"Se\" : "
    "34, \"Br\" : 35, \"Kr\" : 36, \"Rb\" : 37, \"Sr\" : 38, \"Y\"  : 39, "
    "\"Zr\" : 40, \"Nb\" : 41, \"Mo\" : 42, \"Tc\" : 43, \"Ru\" : 44, \"Rh\" : "
    "45, \"Pd\" : 46, \"Ag\" : 47, \"Cd\" : 48, \"In\" : 49, \"Sn\" : 50, "
    "\"Sb\" : 51, \"Te\" : 52, \"I\"  : 53, \"Xe\" : 54, \"Cs\" : 55, \"Ba\" : "
    "56, \"La\" : 57, \"Ce\" : 58, \"Pr\" : 59, \"Nd\" : 60, \"Pm\" : 61, "
    "\"Sm\" : 62, \"Eu\" : 63, \"Gd\" : 64, \"Tb\" : 65, \"Dy\" : 66, \"Ho\" : "
    "67, \"Er\" : 68, \"Tm\" : 69, \"Yb\" : 70, \"Lu\" : 71, \"Hf\" : 72, "
    "\"Ta\" : 73, \"W\"  : 74, \"Re\" : 75, \"Os\" : 76, \"Ir\" : 77, \"Pt\" : "
    "78, \"Au\" : 79, \"Hg\" : 80, \"Tl\" : 81, \"Pb\" : 82, \"Bi\" : 83, "
    "\"Po\" : 84, \"At\" : 85, \"Rn\" : 86, \"Fr\" : 87, \"Ra\" : 88, \"Ac\" : "
    "89, \"Th\" : 90, \"Pa\" : 91, \"U\"  : 92, \"Np\" : 93, \"Pu\" : 94, "
    "\"Am\" : 95, \"Cm\" : 96, \"Bk\" : 97, \"Cf\" : 98, \"Es\" : 99, \"Fm\" : "
    "100, \"Md\" : 101, \"No\" : 102, \"Lr\" : 103, \"Rf\" : 104, \"Ha\" : "
    "105, \"Sg\" : 106, \"Bh\" : 107, \"Hs\" : 108, \"Mt\" : 109 }");

/**************************************************
 *                                                *
 *                parse_lat_to_unita              *
 *                                                *
 **************************************************/

static std::vector<double> parse_lat_to_unita(std::vector<double> lat) {
  double conv = 0.017453292519943; // conv=pi/80.0
  std::vector<double> cang = {lat[3] * conv, lat[4] * conv, lat[5] * conv};
  std::vector<double> gmat = {lat[0] * lat[0],
                              lat[1] * lat[0] * cos(cang[2]),
                              lat[2] * lat[0] * cos(cang[1]),
                              lat[0] * lat[1] * cos(cang[2]),
                              lat[1] * lat[1],
                              lat[2] * lat[1] * cos(cang[0]),
                              lat[0] * lat[2] * cos(cang[1]),
                              lat[1] * lat[2] * cos(cang[0]),
                              lat[2] * lat[2]};
  double deter3 = gmat[0] * (gmat[4] * gmat[8] - gmat[7] * gmat[5]) -
                  gmat[3] * (gmat[1] * gmat[8] - gmat[7] * gmat[2]) +
                  gmat[6] * (gmat[1] * gmat[5] - gmat[4] * gmat[2]);
  double vol = sqrt(deter3);

  double c1 = cos(cang[0]);
  double c2 = cos(cang[1]);
  double s2 = sin(cang[1]);
  double c3 = cos(cang[2]);

  return {lat[0] * s2,
          0.0,
          lat[0] * c2,
          lat[1] * (c3 - c1 * c2) / s2,
          (vol / (lat[0] * lat[2] * s2)),
          lat[1] * c1,
          0.0,
          0.0,
          lat[2]};
}

/**************************************************
 *                                                *
 *                parse_geometry                  *
 *                                                *
 **************************************************/

static json parse_geometry(json geom, int *curptr,
                           std::vector<std::string> lines) {
  json geomjson;

  int cur = *curptr;
  int center = 1;
  int autoz = 0;
  int autosym = 0;
  // double angs_to_au = 1.0/0.52917715;
  double angs_to_au = 1.88972598858;
  double conv = angs_to_au;
  std::vector<std::string> ss;
  std::string geometry = "geometry";

  ss = mystring_split0(lines[cur]);
  if (ss.size() > 1) {
    int nogeomname = mystring_contains(mystring_lowercase(ss[1]), "au") ||
                     mystring_contains(mystring_lowercase(ss[1]), "a.u.") ||
                     mystring_contains(mystring_lowercase(ss[1]), "an") ||
                     mystring_contains(mystring_lowercase(ss[1]), "nm") ||
                     mystring_contains(mystring_lowercase(ss[1]), "na") ||
                     mystring_contains(mystring_lowercase(ss[1]), "pm") ||
                     mystring_contains(mystring_lowercase(ss[1]), "pi") ||
                     mystring_contains(mystring_lowercase(ss[1]), "units") ||
                     mystring_contains(mystring_lowercase(ss[1]), "center") ||
                     mystring_contains(mystring_lowercase(ss[1]), "autoz") ||
                     mystring_contains(mystring_lowercase(ss[1]), "autosym");
    if (!nogeomname)
      geometry = ss[1];
  }

  if (mystring_contains(mystring_lowercase(lines[cur]), " au"))
    conv = 1.0;
  if (mystring_contains(mystring_lowercase(lines[cur]), " a.u."))
    conv = 1.0;
  if (mystring_contains(mystring_lowercase(lines[cur]), " bo"))
    conv = 1.0;
  if (mystring_contains(mystring_lowercase(lines[cur]), " an"))
    conv = angs_to_au;
  if (mystring_contains(mystring_lowercase(lines[cur]), " nm"))
    conv = 10.0 * angs_to_au;
  if (mystring_contains(mystring_lowercase(lines[cur]), " na"))
    conv = 10.0 * angs_to_au;
  if (mystring_contains(mystring_lowercase(lines[cur]), " pm"))
    conv = 0.01 * angs_to_au;
  if (mystring_contains(mystring_lowercase(lines[cur]), " pi"))
    conv = 0.01 * angs_to_au;

  if (mystring_contains(mystring_lowercase(lines[cur]), "nocenter"))
     center = 0;
  else if (mystring_contains(mystring_lowercase(lines[cur]), "center"))
     center = 1;

  if (mystring_contains(mystring_lowercase(lines[cur]), "noautoz"))
     autoz = 0;
  else if (mystring_contains(mystring_lowercase(lines[cur]), "autoz"))
     autoz = 1;

  if (mystring_contains(mystring_lowercase(lines[cur]), "noautosym"))
     autosym = 0;
  else if (mystring_contains(mystring_lowercase(lines[cur]), "autosym"))
     autosym = 1;

  if (mystring_contains(mystring_lowercase(lines[cur]),"symmetry_tolerance"))
     geomjson["symmetry_tolerance"] = std::stod(mystring_split(mystring_lowercase(lines[cur]),"symmetry_tolerance")[1]);

  if (mystring_contains(mystring_lowercase(lines[cur]),"sym_tolerance"))
     geomjson["symmetry_tolerance"] = std::stod(mystring_split(mystring_lowercase(lines[cur]),"sym_tolerance")[1]);


  geomjson["conv"] = conv;
  geomjson["center"] = center;
  geomjson["autoz"] = autoz;
  geomjson["autosym"] = autosym;

  bool is_crystal = false;
  bool is_surface = false;
  bool fractional = false;
  bool twodfractional = false;
  std::vector<double> unita = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  int endcount = 1;
  std::vector<std::string> symbols;
  std::vector<double> coords;
  std::vector<double> velocities;
  std::vector<double> masses;
  std::vector<double> charges;
  ++cur;

  int nion = 0;
  double mm;
  std::string line;

  while (endcount > 0) {
    line = lines[cur];

    if (mystring_contains(line, "system crystal")) {
      is_crystal=true;
      if (mystring_contains(line, "cartesian")) {
        fractional = false;
      } else {
        fractional = true;
      }
      ++endcount;
    } else if (mystring_contains(line, "system surface")) {
      is_surface=true;
      if (mystring_contains(line, "cartesian")) {
        twodfractional = false;
      } else {
        twodfractional = true;
      }
      ++endcount;
    } else if (endcount > 1) {
      if (mystring_contains(line, "lat_a") ||
          mystring_contains(line, "lat_b") ||
          mystring_contains(line, "lat_c") ||
          mystring_contains(line, "alpha") || mystring_contains(line, "beta") ||
          mystring_contains(line, "gamma")) {
        std::vector<double> lat = {0.0, 0.0, 0.0, 90.0, 90.0, 90.0};

        int endsystem = 1;
        while (endsystem > 0) {
          line = mystring_lowercase(lines[cur]);
          ss = mystring_split0(line);
          if (mystring_contains(line, "lat_a"))
            lat[0] = std::stod(ss[1]) * conv;
          if (mystring_contains(line, "lat_b"))
            lat[1] = std::stod(ss[1]) * conv;
          if (mystring_contains(line, "lat_c"))
            lat[2] = std::stod(ss[1]) * conv;
          if (mystring_contains(line, "alpha"))
            lat[3] = std::stod(ss[1]);
          if (mystring_contains(line, "beta"))
            lat[4] = std::stod(ss[1]);
          if (mystring_contains(line, "gamma"))
            lat[5] = std::stod(ss[1]);

          ++cur;
          if (mystring_contains(mystring_lowercase(lines[cur]), "end")) {
            --endsystem;
            --endcount;
          }
        }
        unita = parse_lat_to_unita(lat);

      } else if (mystring_contains(line, "lattice_vectors")) {
        ++cur;
        line = mystring_lowercase(lines[cur]);
        ss = mystring_split0(line);
        if (ss.size() > 2) {
          unita[0] = std::stod(ss[0]) * conv;
          unita[1] = std::stod(ss[1]) * conv;
          unita[2] = std::stod(ss[2]) * conv;
        }
        ++cur;
        line = mystring_lowercase(lines[cur]);
        ss = mystring_split0(line);
        if (ss.size() > 2) {
          unita[3] = std::stod(ss[0]) * conv;
          unita[4] = std::stod(ss[1]) * conv;
          unita[5] = std::stod(ss[2]) * conv;
        }
        ++cur;
        line = mystring_lowercase(lines[cur]);
        ss = mystring_split0(line);
        if (ss.size() > 2) {
          unita[6] = std::stod(ss[0]) * conv;
          unita[7] = std::stod(ss[1]) * conv;
          unita[8] = std::stod(ss[2]) * conv;
        }
      } else if (mystring_contains(line, "sc")) {
        ss = mystring_split0(line);
        if (ss.size() > 1) {
          double a = std::stod(ss[1]) * conv;
          unita = {a, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, a};
        }
      } else if (mystring_contains(line, "fcc")) {
        ss = mystring_split0(line);
        if (ss.size() > 1) {
          double a = 0.5 * std::stod(ss[1]) * conv;
          unita = {a, a, 0.0, a, 0.0, a, 0.0, a, a};
        }
      } else if (mystring_contains(line, "bcc")) {
        ss = mystring_split0(line);
        if (ss.size() > 1) {
          double a = 0.5 * std::stod(ss[1]) * conv;
          unita = {-a, a, a, a, -a, a, a, a, -a};
        }
      }

    }

    // read the geometry
    else if (mystring_trim(line).size() > 1) {

      ss = mystring_split0(lines[cur]);
      symbols.push_back(mystring_capitalize(ss[0]));
      double xx = std::stod(ss[1]);
      double yy = std::stod(ss[2]);
      double zz = std::stod(ss[3]);

      double vx = 0.0;
      double vy = 0.0;
      double vz = 0.0;
      if (ss.size() > 6) {
        bool ffail = false;
        try {
          vx = std::stod(ss[4]);
        } catch (std::exception &ia) {
          vx = 0.0;
          ffail = true;
        }
        try {
          vy = std::stod(ss[5]);
        } catch (std::exception &ia) {
          vy = 0.0;
          ffail = true;
        }
        try {
          vz = std::stod(ss[6]);
        } catch (std::exception &ia) {
          vz = 0.0;
          ffail = true;
        }
        if (ffail) {
          vx = 0.0;
          vy = 0.0;
          vz = 0.0;
        }
      }

      if (fractional) {
        double xxx = unita[0] * xx + unita[3] * yy + unita[6] * zz;
        double yyy = unita[1] * xx + unita[4] * yy + unita[7] * zz;
        double zzz = unita[2] * xx + unita[5] * yy + unita[8] * zz;
        xx = xxx;
        yy = yyy;
        zz = zzz;

        xxx = unita[0] * vx + unita[3] * vy + unita[6] * vz;
        yyy = unita[1] * vx + unita[4] * vy + unita[7] * vz;
        zzz = unita[2] * vx + unita[5] * vy + unita[8] * vz;
        vx = xxx;
        vy = yyy;
        vz = zzz;

      } else if (twodfractional) {
        double xxx = unita[0] * xx + unita[3] * yy + zz*conv;
        double yyy = unita[1] * xx + unita[4] * yy + zz*conv;
        double zzz = unita[2] * xx + unita[5] * yy + zz*conv;
        xx = xxx;
        yy = yyy;
        zz = zzz;

        xxx = unita[0] * vx + unita[3] * vy + vz*conv;
        yyy = unita[1] * vx + unita[4] * vy + vz*conv;
        zzz = unita[2] * vx + unita[5] * vy + vz*conv;
        vx = xxx;
        vy = yyy;
        vz = zzz;
      } else {
        xx *= conv;
        yy *= conv;
        zz *= conv;
        vx *= conv;
        vy *= conv;
        vz *= conv;
      }

      coords.push_back(xx);
      coords.push_back(yy);
      coords.push_back(zz);

      velocities.push_back(vx);
      velocities.push_back(vy);
      velocities.push_back(vz);

      mm = periodic_table_mass[mystring_capitalize(ss[0])];
      if (mystring_contains(mystring_lowercase(line), "mass"))
        mm = std::stod(mystring_split0(mystring_split(line, "mass")[1])[0]);
      masses.push_back(mm);

      mm = (double)periodic_table_Z[mystring_capitalize(ss[0])];
      if (mystring_contains(mystring_lowercase(line), "charge"))
        mm = std::stod(mystring_split0(mystring_split(line, "charge")[1])[0]);
      charges.push_back(mm);

      ++nion;
    }

    ++cur;
    if (mystring_contains(mystring_lowercase(lines[cur]), "end")) {
      --endcount;
      if (endcount > 0)
        ++cur;
    }
  }

  if (center) {
    int ii;
    double xcm = 0.0;
    double ycm = 0.0;
    double zcm = 0.0;
    double tmass = 0.0;
    for (ii = 0; ii < nion; ++ii) {
      xcm += coords[3 * ii] * masses[ii];
      ycm += coords[3 * ii + 1] * masses[ii];
      zcm += coords[3 * ii + 2] * masses[ii];
      tmass += masses[ii];
    }
    // xcm /= ((double) nion);
    // ycm /= ((double) nion);
    // zcm /= ((double) nion);
    xcm /= tmass;
    ycm /= tmass;
    zcm /= tmass;
    for (int ii = 0; ii < nion; ++ii) {
      coords[3 * ii] -= xcm;
      coords[3 * ii + 1] -= ycm;
      coords[3 * ii + 2] -= zcm;
    }
  }

  geomjson["symbols"] = symbols;
  geomjson["coords"] = coords;
  geomjson["velocities"] = velocities;
  geomjson["nion"] = nion;
  geomjson["masses"] = masses;
  geomjson["charges"] = charges;
  geomjson["unita"] = unita;
  geomjson["is_crystal"] = is_crystal;
  geomjson["is_suface"] = is_surface;
  geomjson["fractional"] = fractional;
  geomjson["twodfractional"] = twodfractional;

  *curptr = cur;

  geom[geometry] = geomjson;

  return geom;
}
/**************************************************
 *                                                *
 *              parse_pseudopotentials            *
 *                                                *
 **************************************************/
static json parse_pseudopotentials(json pseudopotentials, int *curptr,
                                   std::vector<std::string> lines) {
  int cur = *curptr;
  int endcount = 1;
  ++cur;
  std::string line;
  std::vector<std::string> ss;

  while (endcount > 0) {
    line = lines[cur];
    if (mystring_contains(lines[cur], "library")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() > 2)
        if (mystring_contains(ss[1], "library"))
          pseudopotentials[ss[0]] = ss[2];
    }
    ++cur;
    if (mystring_contains(lines[cur], "end"))
      --endcount;
  }
  *curptr = cur;
  return pseudopotentials;
}

/**************************************************
 *                                                *
 *                parse_simulation_cell           *
 *                                                *
 **************************************************/
static json parse_simulation_cell(json celljson, int *curptr,
                                  std::vector<std::string> lines) {
  int cur = *curptr;
  int endcount = 1;
  std::string line;
  std::vector<std::string> ss;
  double angs_to_au = 1.88972598858;
  double conv = 1.0;

  if (mystring_contains(mystring_lowercase(lines[cur]), " au"))
    conv = 1.0;
  if (mystring_contains(mystring_lowercase(lines[cur]), " a.u."))
    conv = 1.0;
  if (mystring_contains(mystring_lowercase(lines[cur]), " bo"))
    conv = 1.0;
  if (mystring_contains(mystring_lowercase(lines[cur]), " an"))
    conv = angs_to_au;
  if (mystring_contains(mystring_lowercase(lines[cur]), " nm"))
    conv = 10.0 * angs_to_au;
  if (mystring_contains(mystring_lowercase(lines[cur]), " na"))
    conv = 10.0 * angs_to_au;
  if (mystring_contains(mystring_lowercase(lines[cur]), " pm"))
    conv = 0.01 * angs_to_au;
  if (mystring_contains(mystring_lowercase(lines[cur]), " pi"))
    conv = 0.01 * angs_to_au;

  ++cur;

  while (endcount > 0) {
    line = mystring_lowercase(lines[cur]);
    if (mystring_contains(line, "lattice_vectors")) {
      std::vector<double> unita = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      ++cur;
      line = mystring_lowercase(lines[cur]);
      ss = mystring_split0(line);
      if (ss.size() > 2) {
        unita[0] = std::stod(ss[0]) * conv;
        unita[1] = std::stod(ss[1]) * conv;
        unita[2] = std::stod(ss[2]) * conv;
      }
      ++cur;
      line = mystring_lowercase(lines[cur]);
      ss = mystring_split0(line);
      if (ss.size() > 2) {
        unita[3] = std::stod(ss[0]) * conv;
        unita[4] = std::stod(ss[1]) * conv;
        unita[5] = std::stod(ss[2]) * conv;
      }
      ++cur;
      line = mystring_lowercase(lines[cur]);
      ss = mystring_split0(line);
      if (ss.size() > 2) {
        unita[6] = std::stod(ss[0]) * conv;
        unita[7] = std::stod(ss[1]) * conv;
        unita[8] = std::stod(ss[2]) * conv;
      }
      celljson["unita"] = unita;
    } else if (mystring_contains(line, "lattice")) {
      ++cur;
      line = mystring_lowercase(lines[cur]);
      if (mystring_contains(line, "lat_a") ||
          mystring_contains(line, "lat_b") ||
          mystring_contains(line, "lat_c") ||
          mystring_contains(line, "alpha") || mystring_contains(line, "beta") ||
          mystring_contains(line, "gamma")) {
        std::vector<double> lat = {0.0, 0.0, 0.0, 90.0, 90.0, 90.0};

        int endsystem = 1;
        while (endsystem > 0) {
          line = mystring_lowercase(lines[cur]);
          ss = mystring_split0(line);
          if (mystring_contains(line, "lat_a"))
            lat[0] = std::stod(ss[1]) * conv;
          if (mystring_contains(line, "lat_b"))
            lat[1] = std::stod(ss[1]) * conv;
          if (mystring_contains(line, "lat_c"))
            lat[2] = std::stod(ss[1]) * conv;
          if (mystring_contains(line, "alpha"))
            lat[3] = std::stod(ss[1]);
          if (mystring_contains(line, "beta"))
            lat[4] = std::stod(ss[1]);
          if (mystring_contains(line, "gamma"))
            lat[5] = std::stod(ss[1]);

          ++cur;
          if (mystring_contains(mystring_lowercase(lines[cur]), "end"))
            --endsystem;
        }
        celljson["unita"] = parse_lat_to_unita(lat);
      }
    } else if (mystring_contains(line, "ngrid")) {
      ss = mystring_split0(line);
      if (ss.size() > 3) {
        std::vector<int> ngrid = {std::stoi(ss[1]), std::stoi(ss[2]),
                                  std::stoi(ss[3])};
        celljson["ngrid"] = ngrid;
      }
    } else if (mystring_contains(line, "sc")) {
      ss = mystring_split0(line);
      if (ss.size() > 1) {
        double a = std::stod(ss[1]) * conv;
        std::vector<double> unita = {a, 0.0, 0.0, 0.0, a, 0.0, 0.0, 0.0, a};
        celljson["unita"] = unita;
      }
    } else if (mystring_contains(line, "fcc")) {
      ss = mystring_split0(line);
      if (ss.size() > 1) {
        double a = 0.5 * std::stod(ss[1]) * conv;
        std::vector<double> unita = {a, a, 0.0, a, 0.0, a, 0.0, a, a};
        celljson["unita"] = unita;
      }
    } else if (mystring_contains(line, "bcc")) {
      ss = mystring_split0(line);
      if (ss.size() > 1) {
        double a = 0.5 * std::stod(ss[1]) * conv;
        std::vector<double> unita = {-a, a, a, a, -a, a, a, a, -a};
        celljson["unita"] = unita;
      }
    } else if (mystring_contains(line, "box_delta")) {
      ss = mystring_split0(line);
      if (ss.size() > 1) {
        double a = std::stod(ss[1]);
        celljson["box_delta"] = a * conv;
      }
    } else if (mystring_contains(line, "box_orient")) {
      celljson["box_orient"] = true;
    } else if (mystring_contains(line, "box_different_lengths")) {
      celljson["box_different_lengths"] = true;
    } else if (mystring_contains(line, "boundary_conditions")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        celljson["boundary_conditions"] = ss[1];
      else
        celljson["boundary_conditions"] = "periodic";
    }

    ++cur;
    if (mystring_contains(lines[cur], "end"))
      --endcount;
  }

  *curptr = cur;

  return celljson;
}


/**************************************************
 *                                                *
 *        monkhorst_pack_timereversal_prune       *
 *                                                *
 **************************************************/
 /**
 * @brief Prune and update a JSON array of k-vectors.
 *
 * This function prunes and updates a JSON array of k-vectors (KS) to merge
 * similar k-vectors within a specified tolerance and handle time reversal.
 *
 * @param ks A reference to the JSON array containing k-vectors to be pruned.
 *
 * @details
 * - The function iterates through the k-vectors and merges similar k-vectors
 *   based on a specified tolerance.
 * - K-vectors with their 4th component (weight) less than or equal to 0.0 are
 *   considered invalid and are set to 0.0.
 * - The function also ensures time reversal by negating k-vectors with a
 *   negative norm.
 * - The original JSON array 'ks' is replaced with the updated data.
 *
 * @note
 * - The function assumes that 'ks' is a JSON array of JSON objects, where each
 *   object represents a k-vector with components [kx, ky, kz, weight].
 * - The tolerance for merging similar k-vectors is set to 1.0e-9.
 */
static void monkhorst_pack_timereversal_prune(std::vector<std::vector<double>>& ks) 
{
   size_t nks = ks.size();
   size_t nks2 = 0;
   std::vector<double> kvector;

   std::vector<std::vector<double>> updated_ks;

   for (size_t i = 0; i < nks - 1; ++i)
   {
      if (ks[i][3] > 0.0)
      {
         for (size_t j = i + 1; j < nks; ++j)
         {
            if (std::abs(ks[i][0] + ks[j][0]) < 1.0e-9 &&
                std::abs(ks[i][1] + ks[j][1]) < 1.0e-9 &&
                std::abs(ks[i][2] + ks[j][2]) < 1.0e-9) {
                double tmp3 = ks[i][3] + ks[j][3];
                ks[i][3] = tmp3;
                ks[j][3] = 0.0;
            }
         }
      }
   }

   for (size_t i = 0; i < nks; ++i)
   {
      if (std::abs(ks[i][3]) > 1.0e-9)
      {
         kvector = {ks[i][0], ks[i][1], ks[i][2], ks[i][3]};
         updated_ks.push_back(kvector);
         nks2++;
      }
   }

   for (size_t i = 0; i < nks2; ++i)
   {
      double norm = updated_ks[i][0] + updated_ks[i][1] + updated_ks[i][2];
      if (norm < 0.0) {
        updated_ks[i][0] = -updated_ks[i][0];
        updated_ks[i][1] = -updated_ks[i][1];
        updated_ks[i][2] = -updated_ks[i][2];
      }
   }

   // Replace the original brillouin_zone with the updated data
   ks = updated_ks;
}

// Utility: compare k-points with tolerance
static bool kpoint_equal(const std::vector<double>& k1, const std::vector<double>& k2, double tol=1e-8) {
    return (std::abs(k1[0]-k2[0]) < tol && std::abs(k1[1]-k2[1]) < tol && std::abs(k1[2]-k2[2]) < tol);
}

// Utility: apply a 3x3 rotation matrix to a k-point
static std::vector<double> apply_rotation(const std::vector<double>& k, const std::vector<std::vector<double>>& R) {
    std::vector<double> kout(3,0.0);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            kout[i] += R[i][j]*k[j];
    return kout;
}

// Generate symmetry operations for different point groups
static std::vector<std::vector<std::vector<double>>> generate_symmetry_operations(const std::string& group_name) {
    std::vector<std::vector<std::vector<double>>> ops;
    
    // Identity operation (always present)
    ops.push_back({{1,0,0},{0,1,0},{0,0,1}});
    
    if (group_name == "Oh" || group_name == "O") {
        // Cubic symmetry (Oh or O group)
        // 90, 180, 270 deg rotations about x, y, z
        ops.push_back({{1,0,0},{0,0,-1},{0,1,0}});
        ops.push_back({{1,0,0},{0,-1,0},{0,0,-1}});
        ops.push_back({{1,0,0},{0,0,1},{0,-1,0}});
        ops.push_back({{0,1,0},{-1,0,0},{0,0,1}});
        ops.push_back({{-1,0,0},{0,-1,0},{0,0,1}});
        ops.push_back({{0,-1,0},{1,0,0},{0,0,1}});
        ops.push_back({{0,0,1},{0,1,0},{-1,0,0}});
        ops.push_back({{0,0,-1},{0,1,0},{1,0,0}});
        ops.push_back({{0,0,1},{0,-1,0},{1,0,0}});
        ops.push_back({{0,0,-1},{0,-1,0},{-1,0,0}});
        ops.push_back({{0,1,0},{0,0,1},{1,0,0}});
        ops.push_back({{0,-1,0},{0,0,1},{-1,0,0}});
        ops.push_back({{0,1,0},{0,0,-1},{-1,0,0}});
        ops.push_back({{0,-1,0},{0,0,-1},{1,0,0}});
        ops.push_back({{0,0,1},{1,0,0},{0,1,0}});
        ops.push_back({{0,0,-1},{1,0,0},{0,-1,0}});
        ops.push_back({{0,0,1},{-1,0,0},{0,-1,0}});
        ops.push_back({{0,0,-1},{-1,0,0},{0,1,0}});
        ops.push_back({{1,0,0},{0,1,0},{0,0,-1}});
        ops.push_back({{-1,0,0},{0,1,0},{0,0,1}});
        ops.push_back({{1,0,0},{0,-1,0},{0,0,1}});
        ops.push_back({{-1,0,0},{0,-1,0},{0,0,-1}});
        // 180 deg rotations about face diagonals
        ops.push_back({{0,1,0},{1,0,0},{0,0,-1}});
        ops.push_back({{0,-1,0},{-1,0,0},{0,0,-1}});
        ops.push_back({{0,1,0},{-1,0,0},{0,0,-1}});
        ops.push_back({{0,-1,0},{1,0,0},{0,0,-1}});
        ops.push_back({{0,0,1},{1,0,0},{0,-1,0}});
        ops.push_back({{0,0,-1},{1,0,0},{0,1,0}});
        ops.push_back({{0,0,1},{-1,0,0},{0,1,0}});
        ops.push_back({{0,0,-1},{-1,0,0},{0,-1,0}});
        // 120 deg rotations about body diagonals
        ops.push_back({{0,0,1},{1,0,0},{0,1,0}});
        ops.push_back({{0,0,1},{-1,0,0},{0,-1,0}});
        ops.push_back({{0,0,-1},{1,0,0},{0,-1,0}});
        ops.push_back({{0,0,-1},{-1,0,0},{0,1,0}});
        // 180 deg rotations about axes through edge centers
        ops.push_back({{0,1,0},{0,0,1},{1,0,0}});
        ops.push_back({{0,-1,0},{0,0,1},{-1,0,0}});
        ops.push_back({{0,1,0},{0,0,-1},{-1,0,0}});
        ops.push_back({{0,-1,0},{0,0,-1},{1,0,0}});
        
        // For Oh group, add improper rotations (multiply each by -1)
        if (group_name == "Oh") {
            size_t n = ops.size();
            for (size_t i=0; i<n; ++i) {
                std::vector<std::vector<double>> Rinv(3, std::vector<double>(3));
                for (int a=0; a<3; ++a)
                    for (int b=0; b<3; ++b)
                        Rinv[a][b] = -ops[i][a][b];
                ops.push_back(Rinv);
            }
        }
    }
    else if (group_name == "D4h" || group_name == "D4") {
        // Tetragonal symmetry (D4h or D4 group)
        // 90, 180, 270 deg rotations about z-axis
        ops.push_back({{0,1,0},{-1,0,0},{0,0,1}});
        ops.push_back({{-1,0,0},{0,-1,0},{0,0,1}});
        ops.push_back({{0,-1,0},{1,0,0},{0,0,1}});
        // 180 deg rotations about x and y axes
        ops.push_back({{1,0,0},{0,-1,0},{0,0,-1}});
        ops.push_back({{-1,0,0},{0,1,0},{0,0,-1}});
        // 180 deg rotations about diagonal axes
        ops.push_back({{0,1,0},{1,0,0},{0,0,-1}});
        ops.push_back({{0,-1,0},{-1,0,0},{0,0,-1}});
        
        // For D4h group, add improper rotations
        if (group_name == "D4h") {
            size_t n = ops.size();
            for (size_t i=0; i<n; ++i) {
                std::vector<std::vector<double>> Rinv(3, std::vector<double>(3));
                for (int a=0; a<3; ++a)
                    for (int b=0; b<3; ++b)
                        Rinv[a][b] = -ops[i][a][b];
                ops.push_back(Rinv);
            }
        }
    }
    else if (group_name == "C4v" || group_name == "C4") {
        // C4v or C4 symmetry (common for surfaces)
        // 90, 180, 270 deg rotations about z-axis
        ops.push_back({{0,1,0},{-1,0,0},{0,0,1}});
        ops.push_back({{-1,0,0},{0,-1,0},{0,0,1}});
        ops.push_back({{0,-1,0},{1,0,0},{0,0,1}});
        
        // For C4v, add vertical mirror planes
        if (group_name == "C4v") {
            // Mirror planes through x and y axes
            ops.push_back({{1,0,0},{0,-1,0},{0,0,1}});
            ops.push_back({{-1,0,0},{0,1,0},{0,0,1}});
            // Diagonal mirror planes
            ops.push_back({{0,1,0},{1,0,0},{0,0,1}});
            ops.push_back({{0,-1,0},{-1,0,0},{0,0,1}});
        }
    }
    else if (group_name == "C3v" || group_name == "C3") {
        // C3v or C3 symmetry
        // 120, 240 deg rotations about z-axis
        ops.push_back({{-0.5,0.866025,0},{-0.866025,-0.5,0},{0,0,1}});
        ops.push_back({{-0.5,-0.866025,0},{0.866025,-0.5,0},{0,0,1}});
        
        // For C3v, add vertical mirror planes
        if (group_name == "C3v") {
            // Mirror planes at 0, 60, 120 degrees
            ops.push_back({{1,0,0},{0,-1,0},{0,0,1}});
            ops.push_back({{-0.5,0.866025,0},{0.866025,0.5,0},{0,0,1}});
            ops.push_back({{-0.5,-0.866025,0},{-0.866025,0.5,0},{0,0,1}});
        }
    }
    else if (group_name == "C2v" || group_name == "C2") {
        // C2v or C2 symmetry
        // 180 deg rotation about z-axis
        ops.push_back({{-1,0,0},{0,-1,0},{0,0,1}});
        
        // For C2v, add vertical mirror planes
        if (group_name == "C2v") {
            // Mirror planes through x and y axes
            ops.push_back({{1,0,0},{0,-1,0},{0,0,1}});
            ops.push_back({{-1,0,0},{0,1,0},{0,0,1}});
        }
    }
    else if (group_name == "Cs") {
        // Mirror plane symmetry only
        // Mirror plane through x-axis (most common for surfaces)
        ops.push_back({{1,0,0},{0,-1,0},{0,0,1}});
    }
    else if (group_name == "Ci") {
        // Inversion center only
        ops.push_back({{-1,0,0},{0,-1,0},{0,0,-1}});
    }
    // For C1 (no symmetry) or unknown groups, only identity is used
    
    return ops;
}

// Helper functions for symmetry detection
static bool check_c4_symmetry(const std::vector<double>& coords, const std::vector<double>& masses, double tol) {
    // Check for 4-fold rotation axis along z-direction
    int nion = coords.size() / 3;
    if (nion < 2) return false;
    
    // Count atoms that are mapped to each other under 90-degree rotation
    int mapped_count = 0;
    for (int i = 0; i < nion; ++i) {
        for (int j = 0; j < nion; ++j) {
            if (i != j && masses[i] == masses[j]) {
                // Check if atom j is the 90-degree rotation of atom i
                double x1 = coords[3*i], y1 = coords[3*i+1], z1 = coords[3*i+2];
                double x2 = coords[3*j], y2 = coords[3*j+1], z2 = coords[3*j+2];
                
                // 90-degree rotation: (x,y) -> (-y,x)
                if (std::abs(x2 + y1) < tol && std::abs(y2 - x1) < tol && std::abs(z2 - z1) < tol) {
                    mapped_count++;
                }
            }
        }
    }
    
    return mapped_count >= nion / 4;  // At least 1/4 of atoms should be mapped
}

static bool check_c3_symmetry(const std::vector<double>& coords, const std::vector<double>& masses, double tol) {
    // Check for 3-fold rotation axis along z-direction
    int nion = coords.size() / 3;
    if (nion < 3) return false;
    
    int mapped_count = 0;
    for (int i = 0; i < nion; ++i) {
        for (int j = 0; j < nion; ++j) {
            if (i != j && masses[i] == masses[j]) {
                double x1 = coords[3*i], y1 = coords[3*i+1], z1 = coords[3*i+2];
                double x2 = coords[3*j], y2 = coords[3*j+1], z2 = coords[3*j+2];
                
                // 120-degree rotation: (x,y) -> (-0.5*x - 0.866*y, 0.866*x - 0.5*y)
                double x_rot = -0.5*x1 - 0.866025*y1;
                double y_rot = 0.866025*x1 - 0.5*y1;
                
                if (std::abs(x2 - x_rot) < tol && std::abs(y2 - y_rot) < tol && std::abs(z2 - z1) < tol) {
                    mapped_count++;
                }
            }
        }
    }
    
    return mapped_count >= nion / 3;  // At least 1/3 of atoms should be mapped
}

static bool check_c2_symmetry(const std::vector<double>& coords, const std::vector<double>& masses, double tol) {
    // Check for 2-fold rotation axis along z-direction
    int nion = coords.size() / 3;
    if (nion < 2) return false;
    
    int mapped_count = 0;
    for (int i = 0; i < nion; ++i) {
        for (int j = 0; j < nion; ++j) {
            if (i != j && masses[i] == masses[j]) {
                double x1 = coords[3*i], y1 = coords[3*i+1], z1 = coords[3*i+2];
                double x2 = coords[3*j], y2 = coords[3*j+1], z2 = coords[3*j+2];
                
                // 180-degree rotation: (x,y) -> (-x,-y)
                if (std::abs(x2 + x1) < tol && std::abs(y2 + y1) < tol && std::abs(z2 - z1) < tol) {
                    mapped_count++;
                }
            }
        }
    }
    
    return mapped_count >= nion / 2;  // At least 1/2 of atoms should be mapped
}

static bool check_mirror_symmetry(const std::vector<double>& coords, const std::vector<double>& masses, double tol) {
    // Check for mirror plane perpendicular to x-axis
    int nion = coords.size() / 3;
    if (nion < 2) return false;
    
    int mapped_count = 0;
    for (int i = 0; i < nion; ++i) {
        for (int j = 0; j < nion; ++j) {
            if (i != j && masses[i] == masses[j]) {
                double x1 = coords[3*i], y1 = coords[3*i+1], z1 = coords[3*i+2];
                double x2 = coords[3*j], y2 = coords[3*j+1], z2 = coords[3*j+2];
                
                // Mirror reflection: (x,y) -> (-x,y)
                if (std::abs(x2 + x1) < tol && std::abs(y2 - y1) < tol && std::abs(z2 - z1) < tol) {
                    mapped_count++;
                }
            }
        }
    }
    
    return mapped_count >= nion / 2;  // At least 1/2 of atoms should be mapped
}

// Enhanced symmetry detection for surfaces and adsorbates
static std::string detect_advanced_symmetry(const std::vector<double>& unita, 
                                           const std::vector<double>& coords,
                                           const std::vector<double>& masses,
                                           double sym_tolerance = 1e-6) {
    if (unita.empty() || coords.empty() || masses.empty()) {
        return "C1";  // No symmetry if no geometry info
    }
    
    // Extract lattice vectors
    double a1[3] = {unita[0], unita[1], unita[2]};
    double a2[3] = {unita[3], unita[4], unita[5]};
    double a3[3] = {unita[6], unita[7], unita[8]};
    
    // Calculate lattice parameters
    double a = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
    double b = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
    double c = sqrt(a3[0]*a3[0] + a3[1]*a3[1] + a3[2]*a3[2]);
    
    // Calculate angles
    double cos_alpha = (a2[0]*a3[0] + a2[1]*a3[1] + a2[2]*a3[2]) / (b * c);
    double cos_beta = (a1[0]*a3[0] + a1[1]*a3[1] + a1[2]*a3[2]) / (a * c);
    double cos_gamma = (a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2]) / (a * b);
    
    // Check for surface-like structures (a3 much larger than a1, a2)
    if (c > 2.0 * std::max(a, b)) {
        // This looks like a surface - analyze 2D symmetry
        
        // Check for square surface (C4v symmetry)
        if (std::abs(a-b) < sym_tolerance && std::abs(cos_gamma) < sym_tolerance) {
            // Additional check: look for 4-fold rotation axis in atomic positions
            bool has_c4 = check_c4_symmetry(coords, masses, sym_tolerance);
            if (has_c4) return "C4v";
        }
        
        // Check for hexagonal surface (C3v symmetry)
        if (std::abs(cos_gamma + 0.5) < sym_tolerance) {
            // Additional check: look for 3-fold rotation axis
            bool has_c3 = check_c3_symmetry(coords, masses, sym_tolerance);
            if (has_c3) return "C3v";
        }
        
        // Check for rectangular surface (C2v symmetry)
        if (std::abs(cos_gamma) < sym_tolerance) {
            bool has_c2 = check_c2_symmetry(coords, masses, sym_tolerance);
            if (has_c2) return "C2v";
        }
        
        // Check for oblique surface (Cs symmetry)
        bool has_mirror = check_mirror_symmetry(coords, masses, sym_tolerance);
        if (has_mirror) return "Cs";
        
        return "C1";  // No 2D symmetry
    }
    
    // Check for bulk crystal symmetries
    // Check for cubic symmetry
    if (std::abs(a-b) < sym_tolerance && std::abs(b-c) < sym_tolerance &&
        std::abs(cos_alpha) < sym_tolerance && std::abs(cos_beta) < sym_tolerance && std::abs(cos_gamma) < sym_tolerance) {
        return "Oh";
    }
    
    // Check for tetragonal symmetry
    if (std::abs(a-b) < sym_tolerance && std::abs(cos_alpha) < sym_tolerance && 
        std::abs(cos_beta) < sym_tolerance && std::abs(cos_gamma) < sym_tolerance) {
        return "D4h";
    }
    
    // Check for orthorhombic symmetry
    if (std::abs(cos_alpha) < sym_tolerance && std::abs(cos_beta) < sym_tolerance && std::abs(cos_gamma) < sym_tolerance) {
        return "D2h";
    }
    
    // Default to lower symmetry
    return "C1";
}

// Detect crystal system and symmetry from lattice vectors and atomic positions
static std::string detect_crystal_symmetry(const std::vector<double>& unita, 
                                          const std::vector<double>& coords,
                                          const std::vector<double>& masses,
                                          double sym_tolerance = 1e-6) {
    return detect_advanced_symmetry(unita, coords, masses, sym_tolerance);
}

// General k-point reduction using detected symmetry
static void reduce_kpoints_by_symmetry(std::vector<std::vector<double>>& ks, 
                                      const std::vector<double>& unita = {},
                                      const std::vector<double>& coords = {},
                                      const std::vector<double>& masses = {}) {
    
    // Detect symmetry if lattice information is provided
    std::string group_name = "C1";  // Default to no symmetry
    if (!unita.empty() && !coords.empty() && !masses.empty()) {
        group_name = detect_crystal_symmetry(unita, coords, masses);
        std::cout << "   Detected symmetry: " << group_name << std::endl;
    } else {
        std::cout << "   No geometry information available, using C1 symmetry (no reduction)" << std::endl;
    }
    
    // Generate symmetry operations for the detected group
    auto ops = generate_symmetry_operations(group_name);
    
    // Perform IBZ reduction
    std::vector<std::vector<double>> unique_ks;
    std::vector<double> unique_weights;
    std::vector<bool> used(ks.size(), false);
    double tol = 1e-8;
    size_t nks = ks.size();
    
    for (size_t i = 0; i < nks; ++i) {
        if (used[i]) continue;
        std::vector<double> k = {ks[i][0], ks[i][1], ks[i][2]};
        double w = ks[i][3];
        // Find all symmetry-equivalent k-points
        double total_weight = w;
        used[i] = true;
        for (size_t j = i + 1; j < nks; ++j) {
            if (used[j]) continue;
            std::vector<double> k2 = {ks[j][0], ks[j][1], ks[j][2]};
            bool equiv = false;
            for (const auto& R : ops) {
                std::vector<double> ksym = apply_rotation(k, R);
                // Bring ksym into [-0.5,0.5) range
                for (int d=0; d<3; ++d) {
                    while (ksym[d] >= 0.5) ksym[d] -= 1.0;
                    while (ksym[d] < -0.5) ksym[d] += 1.0;
                }
                if (kpoint_equal(ksym, k2, tol)) {
                    equiv = true;
                    break;
                }
            }
            if (equiv) {
                total_weight += ks[j][3];
                used[j] = true;
            }
        }
        // Store only one representative (the first encountered)
        unique_ks.push_back(k);
        unique_weights.push_back(total_weight);
    }
    
    // Rebuild ks
    ks.clear();
    for (size_t i=0; i<unique_ks.size(); ++i) {
        ks.push_back({unique_ks[i][0], unique_ks[i][1], unique_ks[i][2], unique_weights[i]});
    }
}

// Legacy function for backward compatibility
static void reduce_kpoints_by_cubic_symmetry(std::vector<std::vector<double>>& ks) {
    reduce_kpoints_by_symmetry(ks);
}

/**************************************************
 *                                                *
 *             monkhorst_pack_set                 *
 *                                                *
 **************************************************/
 /**
 * @brief Generate Monkhorst-Pack k-vectors and update a JSON array.
 *
 * This function generates Monkhorst-Pack k-vectors for a given mesh size (nx, ny, nz)
 * and updates a JSON array (ks) with the generated k-vectors. It also handles time reversal
 * if needed.
 *
 * @param nx The number of k-vectors along the x-axis.
 * @param ny The number of k-vectors along the y-axis.
 * @param nz The number of k-vectors along the z-axis.
 * @param ks A reference to the JSON array that will store the generated k-vectors.
 *
 * @details
 * - The function calculates the number of k-vectors (num_kvectors) based on the mesh size.
 * - It calculates the weight for each k-vector based on the number of k-vectors.
 * - The k-vectors are generated and added to the JSON array 'ks' with their weights.
 * - If timereverse is true (all inputs nx, ny, nz are positive), it calls the
 *   'monkhorst_pack_timereversal_prune' function to prune and update the k-vectors
 *   to ensure time reversal symmetry.
 *
 * @note
 * - The generated k-vectors are based on the Monkhorst-Pack grid.
 * - The 'ks' JSON array should be initialized before calling this function.
 */
static void monkhorst_pack_set(const int nx, const int ny, const int nz, std::vector<std::vector<double>>& ks,
                               const std::vector<double>& unita = {},
                               const std::vector<double>& coords = {},
                               const std::vector<double>& masses = {})
{
   int nkx = nx;
   int nky = ny;
   int nkz = nz;
   std::vector<double> kvector;


   bool timereverse = true;
   if (nkx<0)
   {
      nkx = -nkx;
      timereverse = false;
   }
   if (nky<0)
   {
      nky = -nky;
      timereverse = false;
   }
   if (nkz<0)
   {
      nkz = -nkz;
      timereverse = false;
   }
 
   int num_kvectors = nkx*nky*nkz;
   double weight = 1.0/static_cast<double>(num_kvectors);
   double xxx = 1.0/(2.0*nkx);
   double yyy = 1.0/(2.0*nky);
   double zzz = 1.0/(2.0*nkz);
   for (auto i3=0; i3<nkz; ++i3)
   for (auto i2=0; i2<nky; ++i2)
   for (auto i1=0; i1<nkx; ++i1)
   {
      double xx = 1.0 + 2*i1 - nkx;
      double yy = 1.0 + 2*i2 - nky;
      double zz = 1.0 + 2*i3 - nkz;
 
      double kx = xx*xxx;
      double ky = yy*yyy;
      double kz = zz*zzz;
 
      kvector = {kx,ky,kz,weight};
      ks.push_back(kvector);
   }

   size_t initial_kpoints = ks.size();
   std::cout << "   Monkhorst-Pack grid: " << nkx << "x" << nky << "x" << nkz 
             << " = " << initial_kpoints << " k-points" << std::endl;

   if (timereverse) {
      size_t before_timereverse = ks.size();
      monkhorst_pack_timereversal_prune(ks);
      size_t after_timereverse = ks.size();
      std::cout << "   Time-reversal reduction: " << before_timereverse 
                << " -> " << after_timereverse << " k-points" << std::endl;
   }
   
   // Add symmetry reduction using detected symmetry
   size_t before_symmetry = ks.size();
   reduce_kpoints_by_symmetry(ks, unita, coords, masses);
   size_t after_symmetry = ks.size();
   std::cout << "   Symmetry reduction: " << before_symmetry 
             << " -> " << after_symmetry << " k-points" << std::endl;
   std::cout << "   Total reduction: " << initial_kpoints 
             << " -> " << after_symmetry << " k-points (factor of " 
             << std::fixed << std::setprecision(1) << (double)initial_kpoints/after_symmetry << "x)" << std::endl;
}


/**************************************************
 *                                                *
 *                parse_brillouin_zone            *
 *                                                *
 **************************************************/
/**
 * @brief Parse Brillouin zone information from a list of lines.
 *
 * This function parses Brillouin zone information from a list of lines and updates
 * a JSON object (brillouinjson) with the parsed data.
 *
 * @param brillouinjson A JSON object that will store the parsed Brillouin zone information.
 * @param curptr A pointer to the current line index to keep track of parsing progress.
 * @param lines A vector of strings representing the lines of input data to be parsed.
 * @return A JSON object containing the parsed Brillouin zone information.
 *
 * @details
 * - The function iterates through the lines of input data, looking for specific keywords.
 * - If it encounters a "kvector" keyword, it extracts the k-vector components and weight
 *   and adds them to the JSON object 'brillouinjson'.
 * - If it encounters a "monkhorst-pack" keyword, it extracts mesh size information and
 *   calls the 'monkhorst_pack_set' function to generate and update k-vectors.
 * - Other keywords like "path" and "max_kpoints_print" can be added as needed.
 * - The parsing process continues until the "end" keyword is encountered.
 *
 * @note
 * - The 'brillouinjson' JSON object should be initialized before calling this function.
 * - This function assumes that 'lines' contains the input data to be parsed.
 * - Additional keywords and parsing logic can be added as needed for specific input formats.
 */
static json parse_brillouin_zone(json brillouinjson, int *curptr, std::vector<std::string> lines, 
                                 const std::vector<double>& unita = {},
                                 const std::vector<double>& coords = {},
                                 const std::vector<double>& masses = {}) 
{
   int cur = *curptr;
   int endcount = 1;
   ++cur;
   std::string line;
   std::vector<std::string> ss;

   while (endcount > 0) 
   {
      line = mystring_lowercase(lines[cur]);
      if (mystring_contains(line, "kvector")) {
         ss = mystring_split0(line);
         if (ss.size() > 3)
         {
            std::vector<double> kvector;
            if (ss.size() > 4) 
               kvector = {std::stod(ss[1]), std::stod(ss[2]), std::stod(ss[3]), std::stod(ss[4])};
            else
               kvector = {std::stod(ss[1]), std::stod(ss[2]), std::stod(ss[3]), -1.0};
         
            if (brillouinjson["kvectors"].is_null())
            {
               std::vector<std::vector<double>> kvectors;
               kvectors.push_back(kvector);
               brillouinjson["kvectors"] = kvectors;
            }
            else
               brillouinjson["kvectors"].push_back(kvector);
         }
      } else if (mystring_contains(line, "monkhorst-pack")) {
         ss = mystring_split0(line);
         int nkx=1,nky=1,nkz=1;
         if (ss.size() > 1) nkx = std::stoi(ss[1]);
         if (ss.size() > 2) nky = std::stoi(ss[2]);
         if (ss.size() > 3) nkz = std::stoi(ss[3]);

         std::vector<std::vector<double>> kvectors;
         // Check if "kvectors" key exists in the JSON object
         if (!brillouinjson["kvectors"].is_null() && brillouinjson["kvectors"].is_array()) 
         {
             // Convert the JSON array to a vector of vectors of doubles
             kvectors = brillouinjson["kvectors"].get<std::vector<std::vector<double>>>();
         }
         // Use geometry information passed from the parsing context
         monkhorst_pack_set(nkx,nky,nkz,kvectors, unita, coords, masses);
         brillouinjson["kvectors"] = kvectors;

      } else if (mystring_contains(line, "path")) {
         //band_path_set(brillouinjson);

      } else if (mystring_contains(line, "zone_name")) {
         ss = mystring_split0(line);
         if (ss.size() > 1) brillouinjson["zone_name"] = ss[1];

      } else if (mystring_contains(line, "max_kpoints_print")) {
         ss = mystring_split0(line);
         if (ss.size() > 1) brillouinjson["max_kpoints_print"] = std::stoi(ss[1]);
      }

      ++cur;
      if (mystring_contains(lines[cur], "end"))
         --endcount;
   }

   *curptr = cur;

   return brillouinjson;
}


/**************************************************
 *                                                *
 *                parse_steepest_descent          *
 *                                                *
 **************************************************/

static json parse_steepest_descent(json sdjson, int *curptr, std::vector<std::string> lines) 
{
   int cur = *curptr;
   int endcount = 1;
   ++cur;
   std::string line;
   std::vector<std::string> ss;
 
   while (endcount > 0) {
      line = mystring_lowercase(lines[cur]);
     
      if (mystring_contains(line, "loop")) {
         std::vector<int> loop = {1,1};
         ss = mystring_split0(line);
         if (ss.size() > 1) loop[0] = std::stoi(ss[1]);
         if (ss.size() > 2) loop[1] = std::stoi(ss[2]);
         sdjson["loop"] = loop;
      } else if (mystring_contains(line, "xc")) {
         sdjson["xc"] = mystring_trim(mystring_split(line, "xc")[1]);
      } else if (mystring_contains(line, "geometry_optimize")) {
         sdjson["geometry_optimize"] = true;
         if (mystring_contains(line, " off"))   sdjson["geometry_optimize"] = false;
         if (mystring_contains(line, " no"))    sdjson["geometry_optimize"] = false;
         if (mystring_contains(line, " false")) sdjson["geometry_optimize"] = false;
         if (mystring_contains(line, " on"))    sdjson["geometry_optimize"] = true;
         if (mystring_contains(line, " yes"))   sdjson["geometry_optimize"] = true;
         if (mystring_contains(line, " true"))  sdjson["geometry_optimize"] = true;
      } else if (mystring_contains(line, "input_wavefunction_filename")) {
         ss = mystring_split0(line);
         if (ss.size() > 1)
            sdjson["input_wavefunction_filename"] = ss[1];
      } else if (mystring_contains(line, "output_wavefunction_filename")) {
         ss = mystring_split0(line);
         if (ss.size() > 1)
            sdjson["output_wavefunction_filename"] = ss[1];
      } else if (mystring_contains(line, "time_step")) {
         ss = mystring_split0(line);
         if (ss.size() > 1)
            sdjson["time_step"] = std::stod(ss[1]);
      } else if (mystring_contains(line, "fake_mass")) {
         ss = mystring_split0(line);
         if (ss.size() > 1)
            sdjson["fake_mass"] = std::stod(ss[1]);
      } else if (mystring_contains(line, "cutoff")) {
         ss = mystring_split0(line);
         if (ss.size() == 2) sdjson["cutoff"] = {std::stod(ss[1]), 2 * std::stod(ss[1])};
         if (ss.size() > 2)  sdjson["cutoff"] = {std::stod(ss[1]), std::stod(ss[2])};
      } else if (mystring_contains(line, "tolerances")) {
         ss = mystring_split0(line);
         if (ss.size()==2) sdjson["tolerances"] = {std::stod(ss[1]), std::stod(ss[1]), 1.0e-4};
         if (ss.size()==3) sdjson["tolerances"] = {std::stod(ss[1]), std::stod(ss[2]), 1.0e-4};
         if (ss.size()>3)  sdjson["tolerances"] = {std::stod(ss[1]), std::stod(ss[2]), std::stod(ss[3])};
      } else if (mystring_contains(line, "deltae_check")) {
         if (mystring_contains(line, " off"))   sdjson["deltae_check"] = false;
         if (mystring_contains(line, " no"))    sdjson["deltae_check"] = false;
         if (mystring_contains(line, " false")) sdjson["deltae_check"] = false;
         if (mystring_contains(line, " on"))    sdjson["deltae_check"] = true;
         if (mystring_contains(line, " yes"))   sdjson["deltae_check"] = true;
         if (mystring_contains(line, " true"))  sdjson["deltae_check"] = true;
      }
      ++cur;
      if (mystring_contains(lines[cur], "end"))
         --endcount;
   }

   *curptr = cur;

   return sdjson;
}

/**************************************************
 *                                                *
 *                parse_car_parrinello            *
 *                                                *
 **************************************************/

static json parse_car_parrinello(json cpmdjson, int *curptr,
                                 std::vector<std::string> lines) {
  int cur = *curptr;
  int endcount = 1;
  ++cur;
  std::string line;
  std::vector<std::string> ss;

  while (endcount > 0) {
    line = mystring_lowercase(lines[cur]);
    if (mystring_contains(line, "loop")) {
      std::vector<int> loop = {1, 1};
      ss = mystring_split0(line);
      if (ss.size() > 1)
        loop[0] = std::stoi(ss[1]);
      if (ss.size() > 2)
        loop[1] = std::stoi(ss[2]);
      cpmdjson["loop"] = loop;
    } else if (mystring_contains(line, "xc")) {
      cpmdjson["xc"] = mystring_trim(mystring_split(line, "xc")[1]);
    } else if (mystring_contains(line, "input_wavefunction_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["input_wavefunction_filename"] = ss[1];
    } else if (mystring_contains(line, "output_wavefunction_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["output_wavefunction_filename"] = ss[1];
    } else if (mystring_contains(line, "input_v_wavefunction_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["input_v_wavefunction_filename"] = ss[1];
    } else if (mystring_contains(line, "output_v_wavefunction_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["output_v_wavefunction_filename"] = ss[1];
    } else if (mystring_contains(line, "time_step")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["time_step"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "fake_mass")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["fake_mass"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "cutoff")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        cpmdjson["cutoff"] = {std::stod(ss[1]), 2 * std::stod(ss[1])};
      if (ss.size() > 2)
        cpmdjson["cutoff"] = {std::stod(ss[1]), std::stod(ss[2])};
    } else if (mystring_contains(line, "scaling")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        cpmdjson["scaling"] = {std::stod(ss[1]), std::stod(ss[1])};
      if (ss.size() > 2)
        cpmdjson["scaling"] = {std::stod(ss[1]), std::stod(ss[2])};
    } else if (mystring_contains(line, "sa_decay")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        cpmdjson["sa_decay"] = {std::stod(ss[1]), std::stod(ss[1])};
      if (ss.size() > 2)
        cpmdjson["sa_decay"] = {std::stod(ss[1]), std::stod(ss[2])};
    } else if (mystring_contains(line, "xyz_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["xyz_filename"] = ss[1];
    } else if (mystring_contains(line, "emotion_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["emotion_filename"] = ss[1];
    } else if (mystring_contains(line, "ion_motion_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["ion_motion_filename"] = ss[1];
    } else if (mystring_contains(line, "cif_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["cif_filename"] = ss[1];
    } else if (mystring_contains(line, "eigmotion_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["eigmotion_filename"] = ss[1];
    } else if (mystring_contains(line, "omotion_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["omotion_filename"] = ss[1];
    } else if (mystring_contains(line, "hmotion_filename")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["hmotion_filename"] = ss[1];
    } else if (mystring_contains(line, "fei")) {
      cpmdjson["fei_on"] = true;
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["fei_filename"] = ss[1];
    } else if (mystring_contains(line, "dipole_motion")) {
      cpmdjson["dipole_motion_on"] = true;
      ss = mystring_split0(line);
      if (ss.size() > 1)
        cpmdjson["dipole_motion"] = ss[1];
    }

    // initial_velocities temperatue seed
    else if (mystring_contains(line, "intitial_velocities")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        cpmdjson["initial_velocities"] = {std::stod(ss[1]), 12345};
      else if (ss.size() > 2)
        cpmdjson["initial_velocities"] = {std::stod(ss[1]), std::stoi(ss[2])};
      else
        cpmdjson["initial_velocities"] = {298.15, 12345};
    } else if (mystring_contains(line, "com_shift")) {
      if (mystring_contains(line, " off"))
        cpmdjson["com_shift"] = false;
      if (mystring_contains(line, " on"))
        cpmdjson["com_shift"] = true;
    }

    /* Nose-Hoover Pe Te, Pr Tr start/restart nchain mchain */
    else if (mystring_contains(line, "Nose-Hoover")) {
      if (cpmdjson["Nose-Hoover"].is_null()) {
        json apc;
        cpmdjson["Nose-Hoover"] = apc;
      }
      cpmdjson["Nose-Hoover"]["on"] = true;
      double Tr = 298.16;
      double Te = 298.16;
      double Pr = 1200.0;
      double Pe = 1200.0;
      bool nosers = true;
      int mchain = 1;
      int nchain = 1;

      ss = mystring_split0(line);
      int count = 0;
      for (auto i = 1; i < ss.size(); ++i) {
        if (mystring_isfloat(ss[i])) {
          if (count == 0)
            Pe = std::stod(ss[i]);
          if (count == 1)
            Te = std::stod(ss[i]);
          if (count == 2)
            Pr = std::stod(ss[i]);
          if (count == 3)
            Tr = std::stod(ss[i]);
          if (count == 4)
            mchain = std::stoi(ss[i]);
          if (count == 5)
            nchain = std::stoi(ss[i]);
          ++count;
        } else if (mystring_contains(ss[i], "restart"))
          nosers = true;
        else if (mystring_contains(ss[i], "start"))
          nosers = false;
        else if (mystring_contains(ss[i], "clear"))
          nosers = false;
        else if (mystring_contains(ss[i], "initialize"))
          nosers = false;
      }
      if (count == 1) {
        Pr = Pe;
      }
      if (count == 2) {
        Pr = Pe;
        Tr = Te;
      }
      if (count == 3) {
        Tr = Te;
      }

      cpmdjson["Nose-Hoover"]["Te"] = Te;
      cpmdjson["Nose-Hoover"]["Tr"] = Tr;
      cpmdjson["Nose-Hoover"]["Pe"] = Pe;
      cpmdjson["Nose-Hoover"]["Pr"] = Pr;
      cpmdjson["Nose-Hoover"]["nosers"] = nosers;
      cpmdjson["Nose-Hoover"]["Nchain"] = nchain;
      cpmdjson["Nose-Hoover"]["Mchain"] = mchain;
    }

    /* temperature Tion Pion Tion Pelc start/restart nchain mchain */
    else if (mystring_contains(line, "temperature")) {
      if (cpmdjson["Nose-Hoover"].is_null()) {
        json apc;
        cpmdjson["Nose-Hoover"] = apc;
      }
      cpmdjson["Nose-Hoover"]["on"] = true;
      double Tr = 298.16;
      double Te = 298.16;
      double Pr = 1200.0;
      double Pe = 1200.0;
      bool nosers = true;
      int mchain = 1;
      int nchain = 1;

      ss = mystring_split0(line);
      int count = 0;
      for (auto i = 1; i < ss.size(); ++i) {
        if (mystring_isfloat(ss[i])) {
          if (count == 0)
            Tr = std::stod(ss[i]);
          if (count == 1)
            Pr = std::stod(ss[i]);
          if (count == 2)
            Te = std::stod(ss[i]);
          if (count == 3)
            Pe = std::stod(ss[i]);
          if (count == 4)
            nchain = std::stoi(ss[i]);
          if (count == 5)
            mchain = std::stoi(ss[i]);
          ++count;
        } else if (mystring_contains(ss[i], "restart"))
          nosers = true;
        else if (mystring_contains(ss[i], "start"))
          nosers = false;
        else if (mystring_contains(ss[i], "clear"))
          nosers = false;
        else if (mystring_contains(ss[i], "initialize"))
          nosers = false;
      }
      if (count == 1) {
        Te = Tr;
      }
      if (count == 2) {
        Te = Tr;
        Pe = Pr;
      }
      if (count == 3) {
        Pe = Pr;
      }

      cpmdjson["Nose-Hoover"]["Te"] = Te;
      cpmdjson["Nose-Hoover"]["Tr"] = Tr;
      cpmdjson["Nose-Hoover"]["Pe"] = Pe;
      cpmdjson["Nose-Hoover"]["Pr"] = Pr;
      cpmdjson["Nose-Hoover"]["nosers"] = nosers;
      cpmdjson["Nose-Hoover"]["Nchain"] = nchain;
      cpmdjson["Nose-Hoover"]["Mchain"] = mchain;
    }

    /*
          else if (mystring_contains(line,"xc"))
          else if (mystring_contains(line,"tolerances"))
          else if (mystring_contains(line,"geometry_optimize"))
          else if (mystring_contains(line,"energy"))
          else if (mystring_contains(line,"temperature"))
          else if (mystring_contains(line,"initial_velocities"))
          else if (mystring_contains(line,"fei_quench"))
          else if (mystring_contains(line,"sa_decay"))
    */

    ++cur;
    if (mystring_contains(lines[cur], "end"))
      --endcount;
  }

  *curptr = cur;

  return cpmdjson;
}

/**************************************************
 *                                                *
 *                  parse_dplot                   *
 *                                                *
 **************************************************/
static json parse_dplot(json dplot, int *curptr,
                        std::vector<std::string> lines) {
  int cur = *curptr;
  int endcount = 1;
  ++cur;
  std::string line;
  std::vector<std::string> ss;

  while (endcount > 0) {
    line = lines[cur];
    if (mystring_contains(lines[cur], "vectors")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() > 1)
        dplot["vectors"] = ss[1];
    } else if (mystring_contains(lines[cur], "ncell")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() == 4)
        dplot["ncell"] = {std::stoi(ss[1]), std::stoi(ss[3]), std::stoi(ss[3])};
    } else if (mystring_contains(lines[cur], "origin")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() == 4)
        dplot["origin"] = {std::stod(ss[1]), std::stod(ss[3]),
                           std::stod(ss[3])};
    } else if (mystring_contains(lines[cur], "position_tolerance")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() > 1)
        dplot["position_tolerance"] = std::stod(ss[1]);
    } else if (mystring_contains(lines[cur], "orbital2")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() > 3)
        if (mystring_contains(ss[0], "orbital2"))
          dplot["orbital2-" + ss[1] + "-" + ss[2]] = ss[3];
    } else if (mystring_contains(lines[cur], "orbital")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() > 2)
        if (mystring_contains(ss[0], "orbital"))
          dplot["orbital-" + ss[1]] = ss[2];
    } else if (mystring_contains(lines[cur], "density")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() > 2)
        if (mystring_contains(ss[0], "density"))
          dplot["density-" + ss[1]] = ss[2];
    } else if (mystring_contains(lines[cur], "elf")) {
      ss = mystring_split0(lines[cur]);
      if (ss.size() > 2)
        if (mystring_contains(ss[0], "elf"))
          dplot["elf-" + ss[1]] = ss[2];
    } else if (mystring_contains(lines[cur], "limitxyz")) {
      std::vector<std::string> fourlines = {
          mystring_trim(lines[cur]), mystring_trim(lines[cur + 1]),
          mystring_trim(lines[cur + 2]), mystring_trim(lines[cur + 3])};
      dplot["limitxyz"] = fourlines;
      cur += 3;
    }

    ++cur;
    if (mystring_contains(lines[cur], "end"))
      --endcount;
  }
  *curptr = cur;
  return dplot;
}

/**************************************************
 *                                                *
 *                parse_nwpw                      *
 *                                                *
 **************************************************/

static json parse_nwpw(json nwpwjson, int *curptr,
                       std::vector<std::string> lines,
                       const std::vector<double>& unita = {},
                       const std::vector<double>& coords = {},
                       const std::vector<double>& masses = {}) {
  // json nwpwjson;

  int cur = *curptr;
  int endcount = 1;
  ++cur;
  std::string line;
  std::vector<std::string> ss;

  while (endcount > 0) {
    line = mystring_lowercase(lines[cur]);

    /*
          if (mystring_contains(line,"steepest_descent"))
          if (mystring_contains(line,"car-parrinello"))
    */
    if (mystring_contains(line, "simulation_cell")) 
    {
       if (nwpwjson["simulation_cell"].is_null()) {
         json simulation_cell;
         nwpwjson["simulation_cell"] = simulation_cell;
       }
       *curptr = cur;
       nwpwjson["simulation_cell"] = parse_simulation_cell(nwpwjson["simulation_cell"], curptr, lines);
       cur = *curptr;

    } else if (mystring_contains(line, "brillouin_zone")) {
       if (nwpwjson["brillouin_zone"].is_null() || (mystring_contains(line, "reset"))) 
       {
          json brillouin_zone;
          //json kpoints = nlohmann::json::array();
          nwpwjson["brillouin_zone"] = brillouin_zone;
       }
       *curptr = cur;
       nwpwjson["brillouin_zone"] = parse_brillouin_zone(nwpwjson["brillouin_zone"], curptr, lines, unita, coords, masses);
       cur = *curptr;
    } 
    else if (mystring_contains(line, "monkhorst-pack")) 
    {
       if (nwpwjson["brillouin_zone"].is_null())
       {
          //json kpoints = nlohmann::json::array();
          json brillouinjson;
          nwpwjson["brillouin_zone"] = brillouinjson;
       }
       ss = mystring_split0(line);
       int nkx=1,nky=1,nkz=1;
       if (ss.size() > 1) nkx = std::stoi(ss[1]);
       if (ss.size() > 2) nky = std::stoi(ss[2]);
       if (ss.size() > 3) nkz = std::stoi(ss[3]);

       std::vector<std::vector<double>> kvectors;
       if (!nwpwjson["brillouin_zone"]["kvectors"].is_null() && nwpwjson["brillouin_zone"]["kvectors"].is_array()) 
       {
             // Convert the JSON array to a vector of vectors of doubles
             kvectors = nwpwjson["brillouin_zone"]["kvectors"].get<std::vector<std::vector<double>>>();
       }
       // Use the geometry information passed from the calling context
       // If not available, try to get lattice vectors from simulation_cell
       std::vector<double> local_unita = unita;
       if (local_unita.empty() && !nwpwjson["simulation_cell"]["unita"].is_null()) {
           local_unita = nwpwjson["simulation_cell"]["unita"].get<std::vector<double>>();
       }
       
       monkhorst_pack_set(nkx,nky,nkz,kvectors, local_unita, coords, masses);
       nwpwjson["brillouin_zone"]["kvectors"] = kvectors;


    } 
    else if (mystring_contains(line, "pseudopotentials")) 
    {
       if (nwpwjson["pseudopotentials"].is_null()) {
         json pseudopotentials;
         nwpwjson["pseudopotentials"] = pseudopotentials;
       }
       *curptr = cur;
       nwpwjson["pseudopotentials"] = parse_pseudopotentials(nwpwjson["pseudopotentials"], curptr, lines);
       cur = *curptr;
    } else if (mystring_contains(line, "steepest_descent")) {
       if (nwpwjson["steepest_descent"].is_null()) {
         json steepest_descent;
         nwpwjson["steepest_descent"] = steepest_descent;
       }
       *curptr = cur;
       nwpwjson["steepest_descent"] = parse_steepest_descent(nwpwjson["steepest_descent"], curptr, lines);
       cur = *curptr;
    } else if (mystring_contains(line, "car-parrinello")) {
       if (nwpwjson["car-parrinello"].is_null()) {
         json car_parrinello;
         nwpwjson["car-parrinello"] = car_parrinello;
       }
       *curptr = cur;
       nwpwjson["car-parrinello"] =
           parse_car_parrinello(nwpwjson["car-parrinello"], curptr, lines);
       cur = *curptr;
    } else if (mystring_contains(line, "dplot")) {
       if (nwpwjson["dplot"].is_null()) {
         json dplot;
         nwpwjson["dplot"] = dplot;
       }
       *curptr = cur;
       nwpwjson["dplot"] = parse_dplot(nwpwjson["dplot"], curptr, lines);
       cur = *curptr;
    } else if (mystring_contains(line, "initialize_wavefunction")) {
       if (mystring_contains(line, " off"))
         nwpwjson["initialize_wavefunction"] = false;
       else if (mystring_contains(line, " no"))
         nwpwjson["initialize_wavefunction"] = false;
       else if (mystring_contains(line, " false"))
         nwpwjson["initialize_wavefunction"] = false;
      
       else if (mystring_contains(line, " yes"))
         nwpwjson["initialize_wavefunction"] = true;
       else if (mystring_contains(line, " true"))
         nwpwjson["initialize_wavefunction"] = true;
       else if (mystring_contains(line, " on"))
         nwpwjson["initialize_wavefunction"] = true;
       else
         nwpwjson["initialize_wavefunction"] = true;
    } else if (mystring_contains(line, "io_norbs_max")) {
       ss = mystring_split0(line);
       if (ss.size() == 2)
         nwpwjson["io_norbs_max"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "nobalance")) {
       nwpwjson["nobalance"] = true;
    } else if (mystring_contains(line, "nolagrange")) {
       nwpwjson["nolagrange"] = true;
       if (mystring_contains(line, " off"))        nwpwjson["nolagrange"] = false;
       else if (mystring_contains(line, " false")) nwpwjson["nolagrange"] = false;
       else if (mystring_contains(line, " true"))  nwpwjson["nolagrange"] = true;
       else
           nwpwjson["nolagrange"] = true;
    } else if (mystring_contains(line, "use_grid_cmp")) {
       nwpwjson["use_grid_cmp"] = true;
    } else if (mystring_contains(line, "fast_erf")) {
       nwpwjson["fast_erf"] = true;
    } else if (mystring_contains(line, "mapping")) {
       ss = mystring_split0(line);
       if (ss.size() > 1)
          nwpwjson["mapping"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "initial_psi_random_algorithm")) {
       ss = mystring_split0(line);
       if (ss.size() > 1)
          nwpwjson["initial_psi_random_algorithm"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "initial_wavefunction_guess")) {
       ss = mystring_split0(line);
       if (ss.size() > 1)
          nwpwjson["initial_wavefunction_guess"] = ss[1];
    } else if (mystring_contains(line, "tile_factor")) {
       ss = mystring_split0(line);
       if (ss.size() > 1)
          nwpwjson["tile_factor"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "1d-slab")) {
       nwpwjson["mapping"] = 1;
    } else if (mystring_contains(line, "2d-hilbert")) {
       nwpwjson["mapping"] = 2;
    } else if (mystring_contains(line, "2d-hcurve")) {
       nwpwjson["mapping"] = 3;
    } else if (mystring_contains(line, "np_dimensions")) {
       ss = mystring_split0(line);
       if (ss.size() > 2)
         nwpwjson["np_dimensions"] = {std::stoi(ss[1]), std::stoi(ss[2])};
       if (ss.size() > 3)
         nwpwjson["np_dimensions"] = {std::stoi(ss[1]), std::stoi(ss[2]), std::stoi(ss[3])};

    } else if (mystring_contains(line, "loop")) {
       std::vector<int> loop;
       loop.push_back(1);
       loop.push_back(1);
       ss = mystring_split0(line);
       if (ss.size() > 1)
         loop[0] = std::stoi(ss[1]);
       if (ss.size() > 2)
         loop[1] = std::stoi(ss[2]);
       nwpwjson["loop"] = loop;
    } else if (mystring_contains(line, "bo_steps")) {
       std::vector<int> loop;
       loop.push_back(1);
       loop.push_back(1);
       ss = mystring_split0(line);
       if (ss.size() > 1)
         loop[0] = std::stoi(ss[1]);
       if (ss.size() > 2)
         loop[1] = std::stoi(ss[2]);
       nwpwjson["bo_steps"] = loop;
    } else if (mystring_contains(line, "bo_time_step")) {
       ss = mystring_split0(line);
       if (ss.size() > 1)
          nwpwjson["bo_time_step"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "bo_algorithm")) {
       if (mystring_contains(line, " leap-frog"))
          nwpwjson["bo_algorithm"] = 2;
       else if (mystring_contains(line, " velocity-verlet"))
          nwpwjson["bo_algorithm"] = 1;
       else
          nwpwjson["bo_algorithm"] = 0;
    } else if (mystring_contains(line, "xc")) {
        nwpwjson["xc"] = mystring_trim(mystring_split(line, "xc")[1]);
    } else if (mystring_contains(line, "cutoff")) {
       ss = mystring_split0(line);
       if (ss.size() == 2)
          nwpwjson["cutoff"] = {std::stod(ss[1]), 2 * std::stod(ss[1])};
       if (ss.size() > 2)
          nwpwjson["cutoff"] = {std::stod(ss[1]), std::stod(ss[2])};
    } else if (mystring_contains(line, "ewald_ncut")) {
       ss = mystring_split0(line);
       if (ss.size() == 2)
          nwpwjson["ewald_ncut"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "ewald_rcut")) {
       ss = mystring_split0(line);
       if (ss.size() == 2)
          nwpwjson["ewald_rcut"] = std::stod(ss[1]);
    } else if ((mystring_contains(line, "unrestricted")) || (mystring_contains(line, "odft"))) {
       nwpwjson["ispin"] = 2;
    } else if (mystring_contains(line, "restricted")) {
       nwpwjson["ispin"] = 1;
    } else if (mystring_contains(line, "mult")) {
       ss = mystring_split0(line);
       if (ss.size() == 2)
       {
          int mult = std::stoi(ss[1]);
          nwpwjson["mult"] = mult;
          if (mult > 1) nwpwjson["ispin"] = 2;
       }
    } else if (mystring_contains(line, "eprecondition")) {
       if (ss.size() == 2)
          nwpwjson["eprecondition"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "sprecondition")) {
       if (ss.size() == 2)
          nwpwjson["sprecondition"] = std::stod(ss[1]);

    } else if (mystring_contains(line, "tolerances")) {
       ss = mystring_split0(line);
       if (ss.size() == 2)
          nwpwjson["tolerances"] = {std::stod(ss[1]), std::stod(ss[1]), 1.0e-4};
       if (ss.size() == 3)
          nwpwjson["tolerances"] = {std::stod(ss[1]), std::stod(ss[2]), 1.0e-4};
       if (ss.size() > 3)
          nwpwjson["tolerances"] = {std::stod(ss[1]), std::stod(ss[2]),
                                   std::stod(ss[3])};
    } else if (mystring_contains(line, "time_step")) {
       ss = mystring_split0(line);
       if (ss.size() > 1)
          nwpwjson["time_step"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "intitial_velocities")) {
       ss = mystring_split0(line);
       if (ss.size() == 2)
         nwpwjson["initial_velocities"] = {std::stod(ss[1]), 12345};
       else if (ss.size() > 2)
         nwpwjson["initial_velocities"] = {std::stod(ss[1]), std::stoi(ss[2])};
       else
         nwpwjson["initial_velocities"] = {298.15, 12345};


    } else if (mystring_contains(line, "scf")) {
       if (mystring_contains(line, "ks-grassmann-cg"))
          nwpwjson["minimizer"] = 3;
       else if (mystring_contains(line, "ks-grassmann-lmbfgs"))
          nwpwjson["minimizer"] = 6;
       else if (mystring_contains(line, "potential"))
          nwpwjson["minimizer"] = 5;
       else
          nwpwjson["minimizer"] = 8;

       if (mystring_contains(line, "ks-block-cg"))        nwpwjson["ks_algorithm"] = -1;
       if (mystring_contains(line, "ks-cg"))              nwpwjson["ks_algorithm"] = 0;
       if (mystring_contains(line, "ks-rmm-diis"))        nwpwjson["ks_algorithm"] = 1;
       if (mystring_contains(line, "ks-grassmann-cg"))     nwpwjson["ks_algorithm"] = 2;
       if (mystring_contains(line, "ks-grassmann-lmbfgs")) nwpwjson["ks_algorithm"] = 3;

       if (mystring_contains(line, "simple"))  nwpwjson["scf_algorithm"] = 0;
       if (mystring_contains(line, "broyden")) nwpwjson["scf_algorithm"] = 1;
       if (mystring_contains(line, "pulay"))         nwpwjson["scf_algorithm"] = 2;
       if (mystring_contains(line, "johnson-pulay")) nwpwjson["scf_algorithm"] = 2;
       if (mystring_contains(line, "diis"))          nwpwjson["scf_algorithm"] = 2;
       if (mystring_contains(line, "anderson"))      nwpwjson["scf_algorithm"] = 3;
       if (mystring_contains(line, "thomas-fermi"))  nwpwjson["scf_algorithm"] = 4;
       if (mystring_contains(line, "alpha")) 
          nwpwjson["scf_alpha"] = std::stod(mystring_trim(mystring_split(line, "alpha")[1]));
       if (mystring_contains(line, "beta")) 
          nwpwjson["scf_beta"] = std::stod(mystring_trim(mystring_split(line, "beta")[1]));
       if (mystring_contains(line, "kerker")) 
          nwpwjson["kerker_g0"] = std::stod(mystring_trim(mystring_split(line, "kerker")[1]));
       if (mystring_contains(line, " iterations")) 
          nwpwjson["ks_maxit_orb"] = std::stoi(mystring_trim(mystring_split(line, " iterations")[1]));
       if (mystring_contains(line, "outer_iterations")) 
          nwpwjson["ks_maxit_orbs"] = std::stoi(mystring_trim(mystring_split(line, "outer_iterations")[1]));
       if (mystring_contains(line, "diis_histories")) 
          nwpwjson["diis_histories"] = std::stoi(mystring_trim(mystring_split(line, "diis_histories")[1]));
       if (mystring_contains(line, "extra_rotate")) 
       {
          nwpwjson["scf_extra_rotate"] = true;
          if (mystring_contains(line, "off"))   nwpwjson["scf_extra_rotate"] = false;
          if (mystring_contains(line, "no"))    nwpwjson["scf_extra_rotate"] = false;
          if (mystring_contains(line, "false")) nwpwjson["scf_extra_rotate"] = false;
          if (mystring_contains(line, "on"))    nwpwjson["scf_extra_rotate"] = true;
          if (mystring_contains(line, "yes"))   nwpwjson["scf_extra_rotate"] = true;
          if (mystring_contains(line, "true"))  nwpwjson["scf_extra_rotate"] = true;
       }


    // smear 
    //SMEAR <sigma default 0.001> 
    //[TEMPERATURE <temperature>] 
    //[FERMI || GAUSSIAN || MARZARI-VANDERBILT default FERMI] 
    //[ORBITALS <integer orbitals default 4>] 

    } else if (mystring_contains(line, "cg")) {
       if (mystring_contains(line, "stiefel"))
         nwpwjson["minimizer"] = 4;
       else if (mystring_contains(line, "stich"))
         nwpwjson["minimizer"] = 9;
       else
         nwpwjson["minimizer"] = 1;

    } else if (mystring_contains(line, "lmbfgs")) {
       if (mystring_contains(line, "stiefel"))
         nwpwjson["minimizer"] = 7;
       else if (mystring_contains(line, "stich"))
         nwpwjson["minimizer"] = 10;
       else
         nwpwjson["minimizer"] = 2;

       int lmbfgs_size = 2;
       ss = mystring_split0(line);
       for (auto iis = 0; iis < ss.size(); ++iis)
         if (mystring_isfloat(ss[iis]))
            lmbfgs_size = std::stoi(ss[iis]);
       if (lmbfgs_size > 2)
          nwpwjson["lmbfgs_size"] = lmbfgs_size;

    } else if (mystring_contains(line, "smear")) {
       double kT=0.001;
       double kb=3.16679e-6;
       double temperature;

       ss = mystring_split0(line);
       if (ss.size()>1)
         if (mystring_isfloat(ss[1]))
           kT = std::stod(ss[1]);

       temperature = kT/kb;

       nwpwjson["fractional"] = true;
       nwpwjson["fractional_frozen"] = false;
       nwpwjson["fractional_orbitals"] = {4,4};
       nwpwjson["fractional_kT"] = kT;
       nwpwjson["fractional_temperature"] = temperature;
       nwpwjson["fractional_smeartype"]   = 2;

       if (mystring_contains(line, "fixed"))              nwpwjson["fractional_smeartype"] = -1;
       if (mystring_contains(line, "step"))               nwpwjson["fractional_smeartype"] = 0;
       if (mystring_contains(line, "fermi"))              nwpwjson["fractional_smeartype"] = 1;
       if (mystring_contains(line, "gaussian"))           nwpwjson["fractional_smeartype"] = 2;
       if (mystring_contains(line, "hermite"))            nwpwjson["fractional_smeartype"] = 3;
       if (mystring_contains(line, "marzari-vanderbilt")) nwpwjson["fractional_smeartype"] = 4;
       if (mystring_contains(line, "methfessel-paxton"))  nwpwjson["fractional_smeartype"] = 5;
       if (mystring_contains(line, "cold"))               nwpwjson["fractional_smeartype"] = 6;
       if (mystring_contains(line, "lorentzian"))         nwpwjson["fractional_smeartype"] = 7;
       if (mystring_contains(line, "no correction"))      nwpwjson["fractional_smeartype"] = 8;
       if (mystring_contains(line, "frozen"))             nwpwjson["fractional_frozen"] = true;
       if (mystring_contains(line, "orbitals"))  
       {
          std::vector<int> norbs;
          norbs.push_back(1);
          norbs.push_back(1);
          //ss = mystring_split0(line);
          ss = mystring_split0(mystring_trim(mystring_split(line, "orbitals")[1]));
          if (ss.size() > 1)
             norbs[0] = std::stoi(ss[1]);
          if (ss.size() > 2)
             norbs[1] = std::stoi(ss[2]);
          nwpwjson["fractional_orbitals"] = norbs;
       }
       if (mystring_contains(line, "filling"))  
       {
          std::string rr = " " + mystring_ireplace(mystring_split(mystring_split(mystring_trim(mystring_split(line, "filling")[1]),"]")[0],"[")[1], ",", " ");
          //std::vector<double> filling = mystring_double_list(mystring_ireplace(mystring_split(mystring_split(mystring_trim(mystring_split(line, "filling")[1]),"]")[0],"[")[1], ",", " "), " ");
          std::vector<double> filling = mystring_double_list(rr,"");
          nwpwjson["fractional_filling"] = filling;
       }

       if (mystring_contains(line, "alpha ")) 
          nwpwjson["fractional_alpha"] = std::stod(mystring_trim(mystring_split(line, "alpha")[1]));

       if (mystring_contains(line, "alpha_min ")) 
          nwpwjson["fractional_alpha_min"] = std::stod(mystring_trim(mystring_split(line, "alpha_min")[1]));

       if (mystring_contains(line, "alpha_max ")) 
          nwpwjson["fractional_alpha_max"] = std::stod(mystring_trim(mystring_split(line, "alpha_max")[1]));

       if (mystring_contains(line, "beta ")) 
          nwpwjson["fractional_beta"] = std::stod(mystring_trim(mystring_split(line, "beta")[1]));

       if (mystring_contains(line, "gamma ")) 
          nwpwjson["fractional_gamma"] = std::stod(mystring_trim(mystring_split(line, "gamma")[1]));

       if (mystring_contains(line, "temperature"))
       {
          temperature = std::stod(mystring_trim(mystring_split(line, "temperature")[1]));
          kT = temperature/kb;
          nwpwjson["fractional_kT"] = kT;
          nwpwjson["fractional_temperature"] = temperature;
       }

    } else if (mystring_contains(line, "vectors")) {
       if (mystring_contains(line, " input"))
         nwpwjson["input_wavefunction_filename"] = mystring_split0(
             mystring_trim(mystring_split(line, " input")[1]))[0];

       if (mystring_contains(line, " output"))
         nwpwjson["output_wavefunction_filename"] = mystring_split0(
             mystring_trim(mystring_split(line, " output")[1]))[0];

       if (mystring_contains(line, "vinput"))
         nwpwjson["input_v_wavefunction_filename"] = mystring_split0(
             mystring_trim(mystring_split(line, "vinput")[1]))[0];

       if (mystring_contains(line, "voutput"))
         nwpwjson["output_v_wavefunction_filename"] = mystring_split0(
             mystring_trim(mystring_split(line, "voutput")[1]))[0];
    } 
    else if (mystring_contains(line, "translation")) 
    {
       if (mystring_contains(line, " off"))        nwpwjson["fix_translation"] = true;
       else if (mystring_contains(line, " no"))    nwpwjson["fix_translation"] = true;
       else if (mystring_contains(line, " false")) nwpwjson["fix_translation"] = true;
      
       else if (mystring_contains(line, " yes")) nwpwjson["fix_translation"] = false;
       else if (mystring_contains(line, " true")) nwpwjson["fix_translation"] = false;
       else if (mystring_contains(line, " on")) nwpwjson["fix_translation"] = false;
       else if (mystring_contains(line, " allow_translation"))
          nwpwjson["fix_translation"] = false;
       else
          nwpwjson["fix_translation"] = true;
    } 
    else if (mystring_contains(line, "rotation")) 
    {
       if (mystring_contains(line, " off"))        nwpwjson["fix_rotation"] = true;
       else if (mystring_contains(line, " no"))    nwpwjson["fix_rotation"] = true;
       else if (mystring_contains(line, " false")) nwpwjson["fix_rotation"] = true;
       else if (mystring_contains(line, " yes"))   nwpwjson["fix_rotation"] = false;
       else if (mystring_contains(line, " true"))  nwpwjson["fix_rotation"] = false;
       else if (mystring_contains(line, " on"))    nwpwjson["fix_rotation"] = false;
       else
         nwpwjson["fix_rotation"] = false;
    } 
    else if (mystring_contains(line, "apc")) 
    {
       if (nwpwjson["apc"].is_null()) {
         json apc;
         nwpwjson["apc"] = apc;
       }
       nwpwjson["apc"]["on"] = true;
       if (mystring_contains(line, " off"))  nwpwjson["apc"]["on"] = false;
       if (mystring_contains(line, "gc"))    nwpwjson["apc"]["Gc"] = mystring_double_list(line, "gc")[0];
       if (mystring_contains(line, "gamma")) nwpwjson["apc"]["gamma"] = mystring_double_list(line, "gamma");
       if (mystring_contains(line, "u"))     nwpwjson["apc"]["u"] = mystring_double_list(line, "u");
       if (mystring_contains(line, "q"))     nwpwjson["apc"]["q"] = mystring_double_list(line, "q");
    } 
    else if (mystring_contains(line, "born")) 
    {
       if (nwpwjson["born"].is_null()) {
         json born;
         nwpwjson["born"] = born;
       }
       nwpwjson["born"]["on"] = true;
       if (mystring_contains(line, " off"))        nwpwjson["born"]["on"] = false;
       if (mystring_contains(line, " norelax"))    nwpwjson["born"]["relax"] = false;
       else if (mystring_contains(line, " relax")) nwpwjson["born"]["relax"] = true;
       if (mystring_contains(line, "dielec"))
          nwpwjson["born"]["dielec"] = std::stod(mystring_split0(mystring_trim(mystring_split(line, " dielec")[1]))[0]);
       if (mystring_contains(line, "rcut"))
          nwpwjson["born"]["rcut"] = std::stod(mystring_split0(mystring_trim(mystring_split(line, " rcut")[1]))[0]);
       if (mystring_contains(line, "bradii"))
          nwpwjson["born"]["bradii"] = mystring_double_list(line, "bradii");
    } 
    else if (mystring_contains(line, "efield")) 
    {
       if (nwpwjson["efield"].is_null()) {
         json efield;
         nwpwjson["efield"] = efield;
       }
       nwpwjson["efield"]["on"] = true;
       if (mystring_contains(line, " off"))      nwpwjson["efield"]["on"] = false;
       if (mystring_contains(line, " periodic")) nwpwjson["efield"]["type"] = 0;
       if (mystring_contains(line, " apc"))      nwpwjson["efield"]["type"] = 1;
       if (mystring_contains(line, " rgrid"))    nwpwjson["efield"]["type"] = 2;
       if (mystring_contains(line, " center"))   nwpwjson["efield"]["center"] = mystring_double_list(line, "center");
       if (mystring_contains(line, " vector"))   nwpwjson["efield"]["vector"] = mystring_double_list(line, "vector");
      
       if ((nwpwjson["efield"]["on"]) && (mystring_contains(line, " on")) &&
           (!mystring_contains(line, " vector")))
         nwpwjson["efield"]["vector"] = mystring_double_list(line, "on");
    } 
    else if (mystring_contains(line, "generalized_poisson")) 
    {
       std::string check = mystring_split0(mystring_trim(mystring_split(line, "generalized_poisson")[1]))[0];
       nwpwjson["generalized_poisson"]["on"] = true;
       if (mystring_contains(check, "off"))   nwpwjson["generalized_poisson"]["on"] = false;
       if (mystring_contains(check, "no"))    nwpwjson["generalized_poisson"]["on"] = false;
       if (mystring_contains(check, "false")) nwpwjson["generalized_poisson"]["on"] = false;
       if (mystring_contains(check, "on"))    nwpwjson["generalized_poisson"]["on"] = true;
       if (mystring_contains(check, "yes"))   nwpwjson["generalized_poisson"]["on"] = true;
       if (mystring_contains(check, "true"))  nwpwjson["generalized_poisson"]["on"] = true;

       if (mystring_contains(line, " relax_dielec"))
       {
          check = mystring_split0(mystring_trim(mystring_split(line, " relax_dielec")[1]))[0];
          nwpwjson["generalized_poisson"]["relax_dielec"] = true; // default
          if (mystring_contains(check, "off"))   nwpwjson["generalized_poisson"]["relax_dielec"] = false;
          if (mystring_contains(check, "no"))    nwpwjson["generalized_poisson"]["relax_dielec"] = false;
          if (mystring_contains(check, "false")) nwpwjson["generalized_poisson"]["relax_dielec"] = false;
          if (mystring_contains(check, "on"))    nwpwjson["generalized_poisson"]["relax_dielec"] = true;
          if (mystring_contains(check, "yes"))   nwpwjson["generalized_poisson"]["relax_dielec"] = true;
          if (mystring_contains(check, "true"))  nwpwjson["generalized_poisson"]["relax_dielec"] = true;
       }
       if (mystring_contains(line, " fix_dielec"))   nwpwjson["generalized_poisson"]["relax_dielec"] = false;
       if (mystring_contains(line, " unfix_dielec")) nwpwjson["generalized_poisson"]["relax_dielec"] = true;

       if (mystring_contains(line, " cube_dielec"))
       {
          check = mystring_split0(mystring_trim(mystring_split(line, " cube_dielec")[1]))[0];
          nwpwjson["generalized_poisson"]["cube_dielec"] = true;
          if (mystring_contains(check, "off"))   nwpwjson["generalized_poisson"]["cube_dielec"] = false;
          if (mystring_contains(check, "no"))    nwpwjson["generalized_poisson"]["cube_dielec"] = false;
          if (mystring_contains(check, "false")) nwpwjson["generalized_poisson"]["cube_dielec"] = false;
          if (mystring_contains(check, "on"))    nwpwjson["generalized_poisson"]["cube_dielec"] = true;
          if (mystring_contains(check, "yes"))   nwpwjson["generalized_poisson"]["cube_dielec"] = true;
          if (mystring_contains(check, "true"))  nwpwjson["generalized_poisson"]["cube_dielec"] = true;
       }

       if (mystring_contains(line, " dielec"))
          nwpwjson["generalized_poisson"]["dielec"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " dielec")[1]))[0]);
       if (mystring_contains(line, " rho0"))
          nwpwjson["generalized_poisson"]["rho0"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " rho0")[1]))[0]);
       if (mystring_contains(line, " beta"))
          nwpwjson["generalized_poisson"]["beta"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " beta")[1]))[0]);
       if (mystring_contains(line, " rhomin"))
          nwpwjson["generalized_poisson"]["rhomin"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " rhomin")[1]))[0]);
       if (mystring_contains(line, " rhomax"))
          nwpwjson["generalized_poisson"]["rhomax"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " rhomax")[1]))[0]);
       if (mystring_contains(line, " rcut_ion"))
          nwpwjson["generalized_poisson"]["rcut_ion"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " rcut_ion")[1]))[0]);
       if (mystring_contains(line, " alpha"))
          nwpwjson["generalized_poisson"]["alpha"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " alpha")[1]))[0]);
       if (mystring_contains(line, " maxit"))
          nwpwjson["generalized_poisson"]["maxit"] 
          = std::stoi(mystring_split0(mystring_trim(mystring_split(line, " maxit")[1]))[0]);
       if (mystring_contains(line, " rmin"))
          nwpwjson["generalized_poisson"]["rmin"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " rmin")[1]))[0]);
       if (mystring_contains(line, " rmax"))
          nwpwjson["generalized_poisson"]["rmax"] 
          = std::stod(mystring_split0(mystring_trim(mystring_split(line, " rmax")[1]))[0]);
       if (mystring_contains(line, " filter"))
       {
          double x = std::stod(mystring_split0(mystring_trim(mystring_split(line, " filter")[1]))[0]);
          if (std::isnan(x))
             x = 1.0;
          nwpwjson["generalized_poisson"]["filter"] = x;
       }

       if (mystring_contains(line, " andreussi"))  nwpwjson["generalized_poisson"]["model"] = 0;
       if (mystring_contains(line, " andreussi2")) nwpwjson["generalized_poisson"]["model"] = 1;
       if (mystring_contains(line, " fattebert"))  nwpwjson["generalized_poisson"]["model"] = 2;
       if (mystring_contains(line, " sphere"))     nwpwjson["generalized_poisson"]["model"] = 3;
    
    } else if (mystring_contains(line, "staged_gpu_fft")) {
       std::string check = mystring_split0(mystring_trim(mystring_split(line, "staged_gpu_fft")[1]))[0];
       nwpwjson["staged_gpu_fft"]["on"] = true;
       if (mystring_contains(check, "off"))   nwpwjson["staged_gpu_fft"]["on"] = false;
       if (mystring_contains(check, "no"))    nwpwjson["staged_gpu_fft"]["on"] = false;
       if (mystring_contains(check, "false")) nwpwjson["staged_gpu_fft"]["on"] = false;
       if (mystring_contains(check, "on"))    nwpwjson["staged_gpu_fft"]["on"] = true;
       if (mystring_contains(check, "yes"))   nwpwjson["staged_gpu_fft"]["on"] = true;
       if (mystring_contains(check, "true"))  nwpwjson["staged_gpu_fft"]["on"] = true;

    } else if (mystring_contains(line, "fft_container_size")) {
       int nffts_size = 0;
       try {
          nffts_size = std::stoi(mystring_split0(mystring_trim(mystring_split(line, "fft_container_size")[1]))[0]);
       } catch (const std::invalid_argument& e) {
          std::cerr << "Invalid argument: " << e.what() << std::endl;
       } catch (const std::out_of_range& e) {
           std::cerr << "Out of range: " << e.what() << std::endl;
       }
       if ((nffts_size<100) || (nffts_size>0)) 
          nwpwjson["fft_container_size"] = nffts_size;
    }
    else if (mystring_contains(line, "virtual")) {
       std::vector<int> nexcited;
       nexcited.push_back(1);
       nexcited.push_back(1);
       ss = mystring_split0(line);
       if (ss.size() > 1)
         nexcited[0] = std::stoi(ss[1]);
       if (ss.size() > 2)
         nexcited[1] = std::stoi(ss[2]);
       nwpwjson["virtual"] = nexcited;
    }


    ++cur;
    if (mystring_contains(lines[cur], "end"))
      --endcount;
  }

  *curptr = cur;

  return nwpwjson;
}

/**************************************************
 *                                                *
 *                parse_driver                    *
 *                                                *
 **************************************************/
static json parse_driver(json driverjson, int *curptr,
                         std::vector<std::string> lines) {
  // json driverjson;

  int cur = *curptr;
  int endcount = 1;
  ++cur;
  std::string line;
  std::vector<std::string> ss;

  while (endcount > 0) {
    line = mystring_lowercase(lines[cur]);

    if (mystring_contains(line, "maxiter")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        driverjson["maxiter"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "loose")) {
      driverjson["gmax"] = 0.00450;
      driverjson["grms"] = 0.00300;
      driverjson["xmax"] = 0.01800;
      driverjson["xrms"] = 0.01200;
    } else if (mystring_contains(line, "default")) {
      driverjson["gmax"] = 0.00045;
      driverjson["grms"] = 0.00030;
      driverjson["xmax"] = 0.00180;
      driverjson["xrms"] = 0.01200;
    } else if (mystring_contains(line, "tight")) {
      driverjson["gmax"] = 0.000015;
      driverjson["grms"] = 0.00001;
      driverjson["xmax"] = 0.00006;
      driverjson["xrms"] = 0.00004;
    } else if (mystring_contains(line, "gmax")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["gmax"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "grms")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["grms"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "xmax")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["xmax"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "xrms")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["xrms"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "eprec")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["eprec"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "trust")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["trust"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "sadstp")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["tsadstp"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "clear")) {
      driverjson["clear"] = true;
    } else if (mystring_contains(line, "redoautoz")) {
      driverjson["redoautoz"] = true;
    } else if (mystring_contains(line, "inhess")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        driverjson["inhess"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "moddir")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        driverjson["moddir"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "vardir")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        driverjson["vardir"] = std::stoi(ss[1]);
    } else if (mystring_contains(line, "nofirstneg")) {
      driverjson["nofirstneg"] = true;
    } else if (mystring_contains(line, "firstneg")) {
      driverjson["firstneg"] = true;
    } else if (mystring_contains(line, "bscale")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["bscale"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "ascale")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["ascale"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "tscale")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["tscale"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "hscale")) {
      ss = mystring_split0(line);
      if (ss.size() == 2)
        driverjson["hscale"] = std::stod(ss[1]);
    } else if (mystring_contains(line, "noxyz")) {
      driverjson["noxyz"] = true;
    } else if (mystring_contains(line, "xyz")) {
      ss = mystring_split0(line);
      driverjson["xyz"] = ss[1];
    } else if (mystring_contains(line, "lmbfgs_size")) {
      ss = mystring_split0(line);
      if (ss.size() > 1)
        driverjson["lmbfgs_size"] = std::stoi(ss[1]);
    }

    ++cur;
    if (mystring_contains(lines[cur], "end"))
      --endcount;
  }

  *curptr = cur;

  return driverjson;
}

/**************************************************
 *                                                *
 *                parse_constraints               *
 *                                                *
 **************************************************/
 //### Adding reaction ###
//# reaction_type    = AB + C --> AC + B
//# reaction_indexes = 1 2 6
//# reaction_gamma   = -2.400
//##### current gamma=-2.400 ######
//constraints
//   clear
//   spring bondings 1.000000 -2.400 1.0 1 2 -1.0 1 6
//   spring bond 1 2 1.000000 1.2 
//end

static json parse_constraints(json constraintsjson, int *curptr, std::vector<std::string> lines) 
{
   int cur = *curptr;
   int endcount = 1;
   ++cur;
   std::string line;
   std::vector<std::string> ss;
 
   while (endcount > 0) {
     line = mystring_lowercase(lines[cur]);
 
     if (mystring_contains(line, "clear")) {
         constraintsjson.clear();
         
     } else if ((mystring_contains(line, "spring")) && 
                (mystring_contains(line, "bondings"))) {
        int bcount = 0;
        if (!constraintsjson["bondings"].is_null())
           bcount = constraintsjson["bondings"].size();
        ss = mystring_split0(mystring_trim(mystring_split(line, " bondings")[1]));
        int    n0       = (ss.size()-2)/3;
        double Kspring0 = std::stod(ss[0]);
        double gamma0   = std::stod(ss[1]);
        std::vector<int> indexes;
        std::vector<double> coef;
        for (auto i=0; i<n0; ++i)
        {
            coef.push_back(std::stod(ss[3*i+2]));
            indexes.push_back(std::stoi(ss[3*i+3]));
            indexes.push_back(std::stoi(ss[3*i+4]));
        }
        constraintsjson["bondings"][bcount]["coef"]     = coef;
        constraintsjson["bondings"][bcount]["indexes"]  = indexes;
        constraintsjson["bondings"][bcount]["Kspring0"] = Kspring0;
        constraintsjson["bondings"][bcount]["gamma0"]   = gamma0;

     } else if ((mystring_contains(line, "spring")) &&  
                (mystring_contains(line, " cbond"))) {   
        int cbcount = 0;
        if (!constraintsjson["cbond"].is_null())
           cbcount = constraintsjson["cbond"].size();
        ss = mystring_split0(mystring_trim(mystring_split(line, " cbond")[1]));
        int i0 = std::stoi(ss[0]);
        int j0 = std::stoi(ss[1]);
        int k0 = std::stoi(ss[2]);
        int l0 = std::stoi(ss[3]);
        double Kspring0 = std::stod(ss[4]);
        double Rij0 = std::stod(ss[5]);
        double Rkl0 = std::stod(ss[6]);
        constraintsjson["cbond"][cbcount]["i0"]       = i0;
        constraintsjson["cbond"][cbcount]["j0"]       = j0;
        constraintsjson["cbond"][cbcount]["k0"]       = k0;
        constraintsjson["cbond"][cbcount]["l0"]       = l0;
        constraintsjson["cbond"][cbcount]["Kspring0"] = Kspring0;
        constraintsjson["cbond"][cbcount]["Rij0"]     = Rij0;
        constraintsjson["cbond"][cbcount]["Rkl0"]     = Rkl0;
     
     } else if ((mystring_contains(line, "spring")) && 
                (mystring_contains(line, " bond"))) {
        int bcount = 0;
        if (!constraintsjson["bond"].is_null())
           bcount = constraintsjson["bond"].size();
        ss = mystring_split0(mystring_trim(mystring_split(line, " bond")[1]));
        int i0 = std::stoi(ss[0]);
        int j0 = std::stoi(ss[1]);
        double Kspring0 = std::stod(ss[2]);
        double R0 = std::stod(ss[3]);
        constraintsjson["bond"][bcount]["i0"]       = i0;
        constraintsjson["bond"][bcount]["j0"]       = j0;
        constraintsjson["bond"][bcount]["Kspring0"] = Kspring0;
        constraintsjson["bond"][bcount]["R0"]       = R0;
     } else if ((mystring_contains(line, "spring")) &&
                (mystring_contains(line, "angle"))) {
        int acount = 0;
        if (!constraintsjson["angle"].is_null())
           acount = constraintsjson["angle"].size();
        ss = mystring_split0(mystring_trim(mystring_split(line, " angle")[1]));
        int i0 = std::stoi(ss[0]);
        int j0 = std::stoi(ss[1]);
        int k0 = std::stoi(ss[2]);
        double Kspring0 = std::stod(ss[3]);
        double Theta0   = std::stod(ss[4]);
        constraintsjson["angle"][acount]["i0"]       = i0;
        constraintsjson["angle"][acount]["j0"]       = j0;
        constraintsjson["angle"][acount]["k0"]       = k0;
        constraintsjson["angle"][acount]["Kspring0"] = Kspring0;
        constraintsjson["angle"][acount]["Theta0"]   = Theta0;
     } else if ((mystring_contains(line, "spring")) &&
                (mystring_contains(line, " dihedral"))) {
        int dcount = 0;
        if (!constraintsjson["dihedral"].is_null())
           dcount = constraintsjson["dihedral"].size();
        ss = mystring_split0(mystring_trim(mystring_split(line, " dihedral")[1]));
        int i0 = std::stoi(ss[0]);
        int j0 = std::stoi(ss[1]);
        int k0 = std::stoi(ss[2]);
        int l0 = std::stoi(ss[3]);
        double Kspring0 = std::stod(ss[4]);
        double phi0     = std::stod(ss[5]);
        constraintsjson["dihedral"][dcount]["i0"]       = i0;
        constraintsjson["dihedral"][dcount]["j0"]       = j0;
        constraintsjson["dihedral"][dcount]["k0"]       = k0;
        constraintsjson["dihedral"][dcount]["l0"]       = l0;
        constraintsjson["dihedral"][dcount]["Kspring0"] = Kspring0;
        constraintsjson["dihedral"][dcount]["phi0"]     = phi0;
     } else if ((mystring_contains(line, "spring")) &&
                (mystring_contains(line, " cihedral"))) {
        int ccount = 0;
        if (!constraintsjson["cihedral"].is_null())
           ccount = constraintsjson["cihedral"].size();
        ss = mystring_split0(mystring_trim(mystring_split(line, " cihedral")[1]));
        int i0 = std::stoi(ss[0]);
        int j0 = std::stoi(ss[1]);
        int k0 = std::stoi(ss[2]);
        int l0 = std::stoi(ss[3]);
        double Kspring0 = std::stod(ss[4]);
        double phi0     = std::stod(ss[5]);
        constraintsjson["cihedral"][ccount]["i0"]       = i0;
        constraintsjson["cihedral"][ccount]["j0"]       = j0;
        constraintsjson["cihedral"][ccount]["k0"]       = k0;
        constraintsjson["cihedral"][ccount]["l0"]       = l0;
        constraintsjson["cihedral"][ccount]["Kspring0"] = Kspring0;
        constraintsjson["cihedral"][ccount]["phi0"]     = phi0;
     }

 
     ++cur;
     if (mystring_contains(lines[cur], "end"))
       --endcount;
   }
   *curptr = cur;
   return constraintsjson;
}




/**************************************************
 *                                                *
 *                parse_rtdbjson                  *
 *                                                *
 **************************************************/

json parse_rtdbjson(json rtdb) {
  std::vector<std::string> lines = rtdb["nwinput_lines"];
  int n = rtdb["nwinput_nlines"];
  int cur = rtdb["nwinput_cur"];

  bool foundtask = false;
  while ((cur < n) && (!foundtask)) {

    if (mystring_contains(mystring_lowercase(lines[cur]), "unset ")) {
       std::vector<std::string> ss = mystring_split0(lines[cur]);
       if (ss.size() > 1)
          auto count_erase = rtdb.erase(ss[1]);
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "set ")) {
       std::vector<std::string> ss = mystring_split0(lines[cur]);
       if (ss.size() > 2)
          rtdb[ss[1]] = ss[2];
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "start")) {
       rtdb["dbname"] = mystring_trim( mystring_split(mystring_split(lines[cur], "start")[1], "\n")[0]);
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "geometry")) {
       rtdb["geometries"] = parse_geometry(rtdb["geometries"], &cur, lines);
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "title")) {
       rtdb["title"] = mystring_trim(mystring_ireplace( mystring_split(mystring_ireplace(lines[cur], "TITLE", "title"), "title")[1], "\"", ""));
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "charge")) {
       rtdb["charge"] = std::stoi(mystring_trim(mystring_split(mystring_split(lines[cur], "charge")[1], "\n")[0]));
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "nwpw")) {
       // Extract geometry information to pass to nwpw parsing
       std::vector<double> unita, coords, masses;
       if (!rtdb["geometries"].is_null()) {
           // Find the active geometry (usually the first one)
           std::string geom_key = "geometry";
           if (rtdb["geometries"][geom_key].is_null()) {
               // Try to find any geometry key that has coordinates
               for (auto it = rtdb["geometries"].begin(); it != rtdb["geometries"].end(); ++it) {
                   if (it.value().is_object() && !it.value()["coords"].is_null()) {
                       geom_key = it.key();
                       break;
                   }
               }
           }
           
           // If still null, try the first non-null geometry
           if (rtdb["geometries"][geom_key].is_null()) {
               for (auto it = rtdb["geometries"].begin(); it != rtdb["geometries"].end(); ++it) {
                   if (!it.value().is_null() && it.value().is_object()) {
                       geom_key = it.key();
                       break;
                   }
               }
           }
           
           // Debug: print available geometry keys
           std::cout << "Available geometry keys: ";
           for (auto it = rtdb["geometries"].begin(); it != rtdb["geometries"].end(); ++it) {
               std::cout << "\"" << it.key() << "\" ";
           }
           std::cout << std::endl;
           
           if (!rtdb["geometries"][geom_key]["unita"].is_null()) {
               unita = rtdb["geometries"][geom_key]["unita"].get<std::vector<double>>();
               std::cout << "Found unita with " << unita.size() << " elements" << std::endl;
           }
           if (!rtdb["geometries"][geom_key]["coords"].is_null()) {
               coords = rtdb["geometries"][geom_key]["coords"].get<std::vector<double>>();
               std::cout << "Found coords with " << coords.size() << " elements" << std::endl;
           }
           if (!rtdb["geometries"][geom_key]["masses"].is_null()) {
               masses = rtdb["geometries"][geom_key]["masses"].get<std::vector<double>>();
               std::cout << "Found masses with " << masses.size() << " elements" << std::endl;
           }
       }
       rtdb["nwpw"] = parse_nwpw(rtdb["nwpw"], &cur, lines, unita, coords, masses);
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "driver")) {
       rtdb["driver"] = parse_driver(rtdb["driver"], &cur, lines);
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "constraints")) {
       rtdb["constraints"] = parse_constraints(rtdb["constraints"], &cur, lines);
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "task")) {
       rtdb["current_task"] = lines[cur];
       foundtask = true;
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "print")) {
       rtdb["print"] = mystring_trim(
           mystring_split(mystring_split(lines[cur], "print")[1], "\n")[0]);
    } else if (mystring_contains(mystring_lowercase(lines[cur]), "redirect_filename")) {
       rtdb["redirect_filename"] = mystring_trim(mystring_split(mystring_split(lines[cur], "redirect_filename")[1], "\n")[0]);
    }

    ++cur;
  }
  rtdb["nwinput_cur"] = cur;
  rtdb["foundtask"] = foundtask;

  return rtdb;
}

/**************************************************
 *                                                *
 *               parse_nwinput                    *
 *                                                *
 **************************************************/

std::string parse_nwinput(std::string nwinput) 
{
   // fetch the permanent_dir and scratch_dir
   std::string permanent_dir = ".";
   std::string scratch_dir = ".";
   std::string psp_library_dir = "";
   if (mystring_contains(mystring_lowercase(nwinput), "permanent_dir")) {
     if (!mystring_contains( mystring_trim(mystring_split( mystring_split(nwinput, "permanent_dir")[0], "\n") .back()), "#"))
       permanent_dir = mystring_rtrim_slash(mystring_trim(mystring_split( mystring_split(nwinput, "permanent_dir")[1], "\n")[0]));
   }
   if (mystring_contains(mystring_lowercase(nwinput), "scratch_dir")) {
     if (!mystring_contains( mystring_trim( mystring_split(mystring_split(nwinput, "scratch_dir")[0], "\n") .back()), "#"))
       scratch_dir = mystring_rtrim_slash(mystring_trim(
           mystring_split(mystring_split(nwinput, "scratch_dir")[1], "\n")[0]));
   }
   if (mystring_contains(mystring_lowercase(nwinput), "psp_library_dir")) {
     if (!mystring_contains( mystring_trim( mystring_split(mystring_split(nwinput, "psp_library_dir")[0], "\n") .back()), "#"))
        psp_library_dir = mystring_rtrim_slash(mystring_trim(mystring_split(mystring_split(nwinput, "psp_library_dir")[1], "\n")[0]));
   }
 
   // fetch the dbname
   std::string dbname = "nwchemex";
   if (mystring_contains(mystring_lowercase(nwinput), "start")) {
      if (!mystring_contains( mystring_trim(mystring_split( mystring_split(nwinput, "start")[0], "\n") .back()), "#"))
         dbname = mystring_trim(mystring_split(mystring_split(nwinput, "start")[1], "\n")[0]);
   }
 
   json rtdb;
   // read a JSON file
   if ((mystring_contains(mystring_lowercase(nwinput), "restart")) &&
       (!mystring_contains(mystring_trim(mystring_split(mystring_split(nwinput,"restart")[0],"\n").back()),"#"))) 
   {
      
      // read a JSON file
      std::string dbname0 = permanent_dir + "/" + dbname + ".json";
      std::ifstream ifile(dbname0);
      ifile >> rtdb;
   } 
   // intialize the rtdb structure
   else 
   {
      json nwpw, geometries, driver, constraints;
      rtdb["nwpw"] = nwpw;
      rtdb["geometries"] = geometries;
      rtdb["driver"] = driver;
      rtdb["constraints"] = constraints;
   }
 
   // set the dbname
   rtdb["dbname"] = dbname;
 
   // set the permanent_dir and scratch_dir
   rtdb["permanent_dir"] = permanent_dir;
   rtdb["scratch_dir"] = scratch_dir;
   rtdb["psp_library_dir"] = psp_library_dir;
 
   // split nwinput into lines
   std::vector<std::string> lines = mystring_split(nwinput, "\n");
 
   // Remove comments
   for (auto i=lines.begin(); i!=lines.end(); ++i)
      *i = mystring_split(*i, "#")[0];
 
   rtdb["nwinput_lines"] = lines;
   rtdb["nwinput_nlines"] = lines.size();
   rtdb["nwinput_cur"] = 0;
 
   rtdb = parse_rtdbjson(rtdb);
 
   return rtdb.dump();
}

/**************************************************
 *                                                *
 *               parse_rtdbstring                 *
 *                                                *
 **************************************************/

std::string parse_rtdbstring(std::string rtdbstring) {
  auto rtdb = json::parse(rtdbstring);

  rtdb = parse_rtdbjson(rtdb);

  return rtdb.dump();
}

/**************************************************
 *                                                *
 *                 parse_task                     *
 *                                                *
 **************************************************/

int parse_task(std::string rtdbstring) {
  auto rtdb = json::parse(rtdbstring);
  int task = 0;
  if (rtdb["foundtask"]) {
     // Look for pspw jobs
     if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "pspw")) {
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "energy"))           task = 1;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "gradient"))         task = 2;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "optimize"))         task = 3;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "freq"))             task = 4;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "steepest_descent")) task = 5;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "car-parrinello"))   task = 6;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "born-oppenheimer")) task = 7;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "dplot"))            task = 8;
     }
     // Look for band jobs
     if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "band")) {
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "energy"))           task = 11;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "gradient"))         task = 12;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "optimize"))         task = 13;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "freq"))             task = 14;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "steepest_descent")) task = 15;
        if (mystring_contains(mystring_lowercase(rtdb["current_task"]), "steepest_descent")) task = 15;
     }
     // Look for file jobs
     if (mystring_contains(mystring_lowercase(rtdb["current_task"]),"file")) { task=9; }
  }

  return task;
}

/**************************************************
 *                                                *
 *            parse_initialize_wvfnc              *
 *                                                *
 **************************************************/

int parse_initialize_wvfnc(std::string rtdbstring, bool start_bool) {
  auto rtdbjson = json::parse(rtdbstring);
  // bool init_wvfnc = true;
  bool init_wvfnc = start_bool;

  if (rtdbjson["nwpw"]["initialize_wavefunction"].is_boolean())
    init_wvfnc = rtdbjson["nwpw"]["initialize_wavefunction"];

  return init_wvfnc;
}

/**************************************************
 *                                                *
 *       parse_initialize_wvfnc_set               *
 *                                                *
 **************************************************/
std::string parse_initialize_wvfnc_set(std::string rtdbstring, bool setvalue) {
  auto rtdbjson = json::parse(rtdbstring);
  if (rtdbjson["nwpw"]["initialize_wavefunction"].is_boolean())
    rtdbjson["nwpw"]["initialize_wavefunction"] = setvalue;

  return rtdbjson.dump();
}

/**************************************************
 *                                                *
 *       parse_input_wavefunction_filename        *
 *                                                *
 **************************************************/
std::string parse_input_wavefunction_filename(std::string rtdbstring) {
  auto rtdb = json::parse(rtdbstring);

  std::string filename = "eric.movecs";
  if (rtdb["dbname"].is_string()) {
    std::string dbname = rtdb["dbname"];
    filename = dbname + ".movecs";
  }
  // read from nwpw block
  if (rtdb["nwpw"]["input_wavefunction_filename"].is_string())
    filename = rtdb["nwpw"]["input_wavefunction_filename"];

  // append permdir
  std::string permdir = "";
  if (rtdb["permanent_dir"].is_string()) {
    permdir = rtdb["permanent_dir"];
    permdir += "/";
  }

  return permdir + filename;
}

/**************************************************
 *                                                *
 *           parse_gen_lowlevel_rtdbstrs          *
 *                                                *
 **************************************************/
std::vector<std::string> parse_gen_lowlevel_rtdbstrs(std::string rtdbstring) {
  auto rtdbjson = json::parse(rtdbstring);

  double ecut0 = 6000.0;
  double wcut0 = 3000.0;
  if (rtdbjson["nwpw"]["cutoff"][0].is_number_float())
    wcut0 = rtdbjson["nwpw"]["cutoff"][0];
  if (rtdbjson["nwpw"]["cutoff"][1].is_number_float())
    ecut0 = rtdbjson["nwpw"]["cutoff"][1];

  double tol0 = 1.0e-7;
  double tol1 = 1.0e-7;
  if (rtdbjson["nwpw"]["tolerances"][0].is_number_float()) tol0 = rtdbjson["nwpw"]["tolerances"][0];
  if (rtdbjson["nwpw"]["tolerances"][1].is_number_float()) tol1 = rtdbjson["nwpw"]["tolerances"][1];

  int steps;
  double dx;
  if (wcut0 < 5.0) {
    steps = 0;
    dx = 1.0;
  } else if (wcut0 <= 20.0) {
    steps = 2;
    dx = 1.0 / 2.0;
  } else if (wcut0 <= 30.0) {
    steps = 3;
    dx = 1.0 / 3.0;
  } else {
    steps = 4;
    dx = 1.0 / 4.0;
  }

  // run an extra step when task > 3
  if (parse_task(rtdbstring) > 3)
    ++steps;

  // setup steps
  std::vector<std::string> gen_rtdbs;

  double f[3] = {1000.0,100.0,1.0};
  for (int i = 1; i < steps; ++i) {
    rtdbjson["nwpw"]["cutoff"][0] = i * dx * wcut0;
    rtdbjson["nwpw"]["cutoff"][1] = i * dx * ecut0;
    rtdbjson["current_task"] = "energy";

    if (i <= 3)
    {
        rtdbjson["nwpw"]["tolerances"][0] = f[i-1] * tol0;
        rtdbjson["nwpw"]["tolerances"][1] = f[i-1] * tol1;
    }
    else
    {
        rtdbjson["nwpw"]["tolerances"][0] = tol0;
        rtdbjson["nwpw"]["tolerances"][1] = tol1;
    }

    gen_rtdbs.push_back(rtdbjson.dump());
  }

  return gen_rtdbs;
}

/**************************************************
 *                                                *
 *                 parse_write                    *
 *                                                *
 **************************************************/

void parse_write(std::string rtdbstring) {
  auto rtdbjson = json::parse(rtdbstring);
  std::string pdir = rtdbjson["permanent_dir"];
  std::string dbname0 = rtdbjson["dbname"];
  std::cout << "writing rtdbjson = " << pdir + "/" + dbname0 + ".json"
            << std::endl;
  std::ofstream ofile(pdir + "/" + dbname0 + ".json");
  ofile << std::setw(4) << rtdbjson << std::endl;
}

// Extract geometry information from parsed JSON for symmetry detection
static void extract_geometry_info(const json& rtdb, 
                                 std::vector<double>& unita,
                                 std::vector<double>& coords,
                                 std::vector<double>& masses) {
    unita.clear();
    coords.clear();
    masses.clear();
    
    // Extract lattice vectors from simulation cell
    if (!rtdb["nwpw"]["simulation_cell"]["unita"].is_null()) {
        unita = rtdb["nwpw"]["simulation_cell"]["unita"].get<std::vector<double>>();
    }
    
    // Extract atomic coordinates and masses from geometry
    if (!rtdb["geometries"].is_null()) {
        // Find the active geometry (usually the first one)
        std::string geom_key = "geometry";
        if (rtdb["geometries"][geom_key].is_null()) {
            // Try to find any geometry key
            for (auto it = rtdb["geometries"].begin(); it != rtdb["geometries"].end(); ++it) {
                if (it.value().is_object() && !it.value()["coords"].is_null()) {
                    geom_key = it.key();
                    break;
                }
            }
        }
        
        if (!rtdb["geometries"][geom_key]["coords"].is_null()) {
            coords = rtdb["geometries"][geom_key]["coords"].get<std::vector<double>>();
        }
        
        if (!rtdb["geometries"][geom_key]["masses"].is_null()) {
            masses = rtdb["geometries"][geom_key]["masses"].get<std::vector<double>>();
        }
    }
}

} // namespace pwdft
