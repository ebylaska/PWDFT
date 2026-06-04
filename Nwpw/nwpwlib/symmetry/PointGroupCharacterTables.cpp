#include "PointGroupCharacterTables.hpp"
#include <utility> // std::move
#include <iostream>
#include <cmath>

namespace pwdft {

PointGroupCharacterTable PointGroupCharacterTable::build(
    std::string group,
    int order,
    std::vector<std::string> class_names,
    std::vector<int> class_sizes,
    std::vector<Irrep> irreps)
{
    PointGroupCharacterTable t;
    t.group_ = std::move(group);
    t.order_ = order;
    t.class_names_ = std::move(class_names);
    t.class_sizes_ = std::move(class_sizes);
    t.irreps_ = std::move(irreps);
    return t;
}

int PointGroupCharacterTable::irrep_index(const std::string& name) const
{
    for (int i = 0; i < (int)irreps_.size(); ++i)
        if (irreps_[i].name == name) return i;
    throw std::runtime_error("PointGroupCharacterTable: unknown irrep: " + name);
}

// ---- small helpers: build tables ----

static PointGroupCharacterTable make_C1()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "C1", 1,
        {"E"}, {1},
        { Ir{"A", 1, {1}} }
    );
}

static PointGroupCharacterTable make_Cs()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "Cs", 2,
        {"E", "s"}, {1, 1},
        { Ir{"A'",  1, { 1,  1}},
          Ir{"A''", 1, { 1, -1}} }
    );
}

static PointGroupCharacterTable make_Ci()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "Ci", 2,
        {"E", "i"}, {1, 1},
        { Ir{"Ag", 1, { 1,  1}},
          Ir{"Au", 1, { 1, -1}} }
    );
}

static PointGroupCharacterTable make_C2()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "C2", 2,
        {"E", "C2"}, {1, 1},
        { Ir{"A", 1, { 1,  1}},
          Ir{"B", 1, { 1, -1}} }
    );
}

static PointGroupCharacterTable make_C2v()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "C2v", 4,
        {"E", "C2", "sv(xz)", "sv(yz)"},
        {  1,   1,      1,        1   },
        { Ir{"A1", 1, { 1,  1,  1,  1}},
          Ir{"A2", 1, { 1,  1, -1, -1}},
          Ir{"B1", 1, { 1, -1,  1, -1}},
          Ir{"B2", 1, { 1, -1, -1,  1}} }
    );
}

static PointGroupCharacterTable make_D2h()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "D2h", 8,
        {"E","C2(z)","C2(y)","C2(x)","i","sh(xy)","sv(xz)","sv(yz)"},
        {  1,   1,      1,      1,   1,    1,      1,       1 },
        { Ir{"Ag",  1, {1, 1, 1, 1, 1, 1, 1, 1}},
          Ir{"B1g", 1, {1, 1,-1,-1, 1, 1,-1,-1}},
          Ir{"B2g", 1, {1,-1, 1,-1, 1,-1, 1,-1}},
          Ir{"B3g", 1, {1,-1,-1, 1, 1,-1,-1, 1}},
          Ir{"Au",  1, {1, 1, 1, 1,-1,-1,-1,-1}},
          Ir{"B1u", 1, {1, 1,-1,-1,-1,-1, 1, 1}},
          Ir{"B2u", 1, {1,-1, 1,-1,-1, 1,-1, 1}},
          Ir{"B3u", 1, {1,-1,-1, 1,-1, 1, 1,-1}} }
    );
}


static PointGroupCharacterTable make_D4h()
{
    using Ir = PointGroupCharacterTable::Irrep;

    return PointGroupCharacterTable::build(
        "D4h", 16,
        {"E","2C4","C2","2C2'","2C2''","i","2S4","sh","2sv","2sd"},
        {  1,   2,    1,   2,     2,    1,   2,   1,    2,    2 },

        {
            // Fill these with the standard D4h characters (10 entries each).
            Ir{"A1g",1,{1,1,1,1,1, 1,1,1,1,1}},
            Ir{"A2g",1,{1,1,1,-1,-1, 1,1,1,-1,-1}},
            Ir{"B1g",1,{1,-1,1,1,-1, 1,-1,1,1,-1}},
            Ir{"B2g",1,{1,-1,1,-1,1, 1,-1,1,-1,1}},

            Ir{"Eg", 2,{2,0,-2,0,0, 2,0,-2,0,0}},

            Ir{"A1u",1,{1,1,1,1,1, -1,-1,-1,-1,-1}},
            Ir{"A2u",1,{1,1,1,-1,-1, -1,-1,-1,1,1}},
            Ir{"B1u",1,{1,-1,1,1,-1, -1,1,-1,-1,1}},
            Ir{"B2u",1,{1,-1,1,-1,1, -1,1,-1,1,-1}},
            Ir{"Eu", 2,{2,0,-2,0,0, -2,0,2,0,0}},
        }
    );
}


static PointGroupCharacterTable make_D6h()
{
    using Ir = PointGroupCharacterTable::Irrep;

    return PointGroupCharacterTable::build(
        "D6h", 24,
        {"E","2C6","2C3","C2","3C2'","3C2''","i","2S3","2S6","sh","3sd","3sv"},
        {  1,   2,    2,   1,    3,      3,    1,   2,    2,    1,    3,    3 },

        {
          Ir{"A1g",1,{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
          Ir{"A2g",1,{ 1, 1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1}},
          Ir{"B1g",1,{ 1,-1, 1,-1, 1,-1, 1,-1, 1, 1, 1,-1}},
          Ir{"B2g",1,{ 1,-1, 1,-1,-1, 1, 1,-1, 1, 1,-1, 1}},
          Ir{"E1g",2,{ 2, 1,-1,-2, 0, 0, 2, 1,-1,-2, 0, 0}}, // ← FIXED
          Ir{"E2g",2,{ 2,-1,-1, 2, 0, 0, 2,-1,-1, 2, 0, 0}},

          Ir{"A1u",1,{ 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1}},
          Ir{"A2u",1,{ 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1}},
          Ir{"B1u",1,{ 1,-1, 1,-1, 1,-1,-1, 1,-1,-1,-1, 1}},
          Ir{"B2u",1,{ 1,-1, 1,-1,-1, 1,-1, 1,-1,-1, 1,-1}},
          Ir{"E1u",2,{ 2, 1,-1,-2, 0, 0,-2,-1, 1,-2, 0, 0}},
          Ir{"E2u",2,{ 2,-1,-1, 2, 0, 0,-2, 1, 1,-2, 0, 0}}
        }
    );
}

static PointGroupCharacterTable make_D8h()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // D8h: n=8, order = 4n = 32
    // Split (13-class) convention.
    constexpr int n = 8;
    constexpr int order = 4 * n;   // 32
    constexpr int m = n / 2;       // 4
    constexpr double pi = 3.14159265358979323846;

    // Class order (13 classes):
    //  0  E            (1)
    //  1  2C8          (2)   k=1
    //  2  2C8^2        (2)   k=2
    //  3  2C8^3        (2)   k=3
    //  4  C8^4         (1)   k=4
    //  5  4C2'         (4)
    //  6  4C2''        (4)
    //  7  i            (1)
    //  8  2S8          (2)   k=1
    //  9  2S8^3        (2)   k=3
    // 10  sh           (1)
    // 11  4sv          (4)
    // 12  4sd          (4)

    std::vector<std::string> classNames = {
        "E","2C8","2C8^2","2C8^3","C8^4","4C2'","4C2''","i","2S8","2S8^3","sh","4sv","4sd"
    };
    std::vector<int> classSizes = {
        1,  2,    2,     2,     1,    4,     4,    1,  2,    2,     1,   4,    4
    };

    // 1D reps: determined by signs on generators:
    //   r = C8 (principal rotation), s = a chosen C2' (in-plane), inv = i
    // For 1D reps when n is even:
    //   r -> +1 for A-types, r -> -1 for B-types
    //   s -> +1 for "...1" types, s -> -1 for "...2" types (A1 vs A2, B1 vs B2)
    // Then:
    //   C8^k: (r)^k
    //   C2' : s
    //   C2'': (r)^(m) * s   because C2'' = r^{n/4} s is conjugate-shifted; for n=8, r^{2}
    //           but in 1D reps conjugation is irrelevant so it becomes multiplication by r^{2}
    //   i: inv
    //   S8^k = i*C8^k: inv*(r)^k
    //   sh = i*C8^4: inv*(r)^4
    //   sv = i*C2': inv*s
    //   sd = i*C2'': inv*(r)^m*s  with m=4? careful:
    //        here C8^4 is the principal C2 (180°). For the shift between C2' and C2'',
    //        the factor is r^{n/4} = r^2, and (r^2) in 1D is (rSign)^2 = +1 always.
    //        So for n=8, C2' and C2'' have the same 1D character.
    //        (They differ in 2D E-type reps, not in 1D.)
    auto make1D = [&](std::string name, int rSign, int sSign, int invSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(13);

        auto rpow = [&](int k)->double {
            if (rSign == 1) return 1.0;
            return (k % 2) ? -1.0 : 1.0; // (-1)^k
        };

        chi.push_back(1.0);           // E
        chi.push_back(rpow(1));       // 2C8
        chi.push_back(rpow(2));       // 2C8^2
        chi.push_back(rpow(3));       // 2C8^3
        chi.push_back(rpow(4));       // C8^4

        chi.push_back(double(sSign)); // 4C2'
        chi.push_back(double(sSign)); // 4C2''  (same for 1D reps in D8h)

        chi.push_back(double(invSign));            // i
        chi.push_back(double(invSign) * rpow(1));  // 2S8
        chi.push_back(double(invSign) * rpow(3));  // 2S8^3
        chi.push_back(double(invSign) * rpow(4));  // sh

        chi.push_back(double(invSign) * double(sSign)); // 4sv
        chi.push_back(double(invSign) * double(sSign)); // 4sd (same for 1D reps)

        return Ir{name, 1, chi};
    };

    // 2D reps: E_l, l=1..(m-1)=3
    // Proper rotations:
    //   χ(C8^k) = 2 cos(2π l k / 8)
    //   χ(C8^4) = 2 cos(π l) = 2(-1)^l
    // In-plane C2 axes:
    //   For E_l, characters on C2' and C2'': 0 when l is odd, and ±2 when l is even?
    // More precisely (standard D_n character theory):
    //   For E_l (2D), χ(C2') = 0 and χ(C2'') = 0 for all l when using the usual real 2D reps
    //   of rotations in the plane. However D_nh tables for n divisible by 4 show nonzero
    //   values for E_{n/2} only, which is not present here (l runs to 3).
    // For D8h with l=1,2,3: χ on 4C2' and 4C2'' are 0.
    //
    // Improper operations:
    //   g/u parity sets χ(i) = ±2
    //   χ(S8^k) = (±) * 2 cos(2π l k / 8)
    //   χ(sh) = (±) * 2(-1)^l
    //   χ(sv)=χ(sd)=0
    auto make2D = [&](std::string name, int l, int invSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(13);

        auto rot = [&](int k)->double {
            return 2.0 * std::cos(2.0 * pi * double(l) * double(k) / double(n));
        };

        chi.push_back(2.0);     // E
        chi.push_back(rot(1));  // 2C8
        chi.push_back(rot(2));  // 2C8^2
        chi.push_back(rot(3));  // 2C8^3
        double c4 = 2.0 * ((l % 2) ? -1.0 : 1.0); // 2(-1)^l
        chi.push_back(c4);      // C8^4

        chi.push_back(0.0);     // 4C2'
        chi.push_back(0.0);     // 4C2''

        chi.push_back(2.0 * double(invSign)); // i
        chi.push_back(double(invSign) * rot(1)); // 2S8
        chi.push_back(double(invSign) * rot(3)); // 2S8^3
        chi.push_back(double(invSign) * c4);     // sh

        chi.push_back(0.0);     // 4sv
        chi.push_back(0.0);     // 4sd

        return Ir{name, 2, chi};
    };

    std::vector<Ir> irreps;
    irreps.reserve(14);

    // 1D g
    irreps.push_back(make1D("A1g", +1, +1, +1));
    irreps.push_back(make1D("A2g", +1, -1, +1));
    irreps.push_back(make1D("B1g", -1, +1, +1));
    irreps.push_back(make1D("B2g", -1, -1, +1));

    // 2D g: E1g..E3g
    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l) + "g", l, +1));

    // 1D u
    irreps.push_back(make1D("A1u", +1, +1, -1));
    irreps.push_back(make1D("A2u", +1, -1, -1));
    irreps.push_back(make1D("B1u", -1, +1, -1));
    irreps.push_back(make1D("B2u", -1, -1, -1));

    // 2D u: E1u..E3u
    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l) + "u", l, -1));

    return PointGroupCharacterTable::build(
        "D8h",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_D10h()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // D10h: n=10, order = 4n = 40
    // Split-class convention appropriate for even n not divisible by 4:
    // classes: E; 2C10^k (k=1..4); C10^5; 10C2'; i; 2S10^k (k=1..4); sh; 10sv
    // total classes = 1 + 4 + 1 + 1 + 1 + 4 + 1 + 1 = 14
    constexpr int n = 10;
    constexpr int order = 4 * n;   // 40
    constexpr int m = n / 2;       // 5
    constexpr double pi = 3.14159265358979323846;

    std::vector<std::string> classNames;
    std::vector<int> classSizes;

    classNames.push_back("E");      classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) { // k=1..4
        classNames.push_back("2C10^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("C10^5");  classSizes.push_back(1);
    classNames.push_back("10C2'");  classSizes.push_back(n);

    classNames.push_back("i");      classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) { // k=1..4
        classNames.push_back("2S10^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("sh");     classSizes.push_back(1);
    classNames.push_back("10sv");   classSizes.push_back(n);

    // 1D reps via generator signs:
    // r = C10, s = C2' (perpendicular 2-fold), inv = i.
    // 1D reps: r -> +1 (A) or -1 (B); s -> +1 ("1") or -1 ("2"); inv -> ±1 (g/u)
    auto make1D = [&](std::string name, int rSign, int sSign, int invSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(14);

        auto rpow = [&](int k)->double {
            if (rSign == 1) return 1.0;
            return (k % 2) ? -1.0 : 1.0; // (-1)^k
        };

        chi.push_back(1.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rpow(k)); // 2C10^k

        double c5 = rpow(5);                 // C10^5
        chi.push_back(c5);

        chi.push_back(double(sSign));        // 10C2'

        chi.push_back(double(invSign));      // i

        for (int k = 1; k <= m - 1; ++k)     // 2S10^k = i*C10^k
            chi.push_back(double(invSign) * rpow(k));

        chi.push_back(double(invSign) * c5); // sh = i*C10^5
        chi.push_back(double(invSign) * double(sSign)); // 10sv = i*C2'

        return Ir{name, 1, chi};
    };

    // 2D reps E_l, l = 1..(m-1)=4:
    // χ(C10^k)=2cos(2π l k / 10), χ(C10^5)=2cos(π l)=2(-1)^l, χ(C2')=0
    // g/u: χ(i)=±2, χ(S10^k)=±χ(C10^k), χ(sh)=±χ(C10^5), χ(sv)=0
    auto make2D = [&](std::string name, int l, int invSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(14);

        auto rot = [&](int k)->double {
            return 2.0 * std::cos(2.0 * pi * double(l) * double(k) / double(n));
        };

        chi.push_back(2.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rot(k)); // 2C10^k

        double c5 = 2.0 * ((l % 2) ? -1.0 : 1.0); // 2(-1)^l
        chi.push_back(c5);                        // C10^5

        chi.push_back(0.0);                       // 10C2'

        chi.push_back(2.0 * double(invSign));     // i

        for (int k = 1; k <= m - 1; ++k)          // 2S10^k
            chi.push_back(double(invSign) * rot(k));

        chi.push_back(double(invSign) * c5);      // sh
        chi.push_back(0.0);                       // 10sv

        return Ir{name, 2, chi};
    };

    std::vector<Ir> irreps;
    irreps.reserve(4 + (m - 1) + 4 + (m - 1)); // 14

    // 1D g
    irreps.push_back(make1D("A1g", +1, +1, +1));
    irreps.push_back(make1D("A2g", +1, -1, +1));
    irreps.push_back(make1D("B1g", -1, +1, +1));
    irreps.push_back(make1D("B2g", -1, -1, +1));

    // 2D g: E1g..E4g
    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l) + "g", l, +1));

    // 1D u
    irreps.push_back(make1D("A1u", +1, +1, -1));
    irreps.push_back(make1D("A2u", +1, -1, -1));
    irreps.push_back(make1D("B1u", -1, +1, -1));
    irreps.push_back(make1D("B2u", -1, -1, -1));

    // 2D u: E1u..E4u
    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l) + "u", l, -1));

    return PointGroupCharacterTable::build(
        "D10h",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_D20h()
{
    using Ir = PointGroupCharacterTable::Irrep;

    constexpr int n = 20;
    constexpr int order = 4 * n;     // 80
    constexpr int m = n / 2;         // 10
    constexpr double pi = 3.14159265358979323846;

    // Class order we use (23 classes):
    // 0: E
    // 1..9:  2C20^k  (k=1..9)
    // 10:    C20^10
    // 11:    20C2'
    // 12:    i
    // 13..21: 2S20^k (k=1..9)
    // 22:    sh
    // 23:    20sv
    // But that's 24 indices; we actually have 23 classes, so we must combine sh and 20sv as 2 classes:
    // Let's define indices exactly as:
    // 0 E
    // 1..9  Ck   (each is "2C20^k")
    // 10    C10  (single)
    // 11    C2'  (20 elements)
    // 12    i
    // 13..21 Sk  (each is "2S20^k")
    // 22    sh
    // 23    sv   (20 elements)
    // That would be 24 classes, which is wrong.
    //
    // Correct counting: D_nh has classes:
    // E (1)
    // 2C_n^k for k=1..(n/2-1) => 2C20^k for k=1..9 (9)
    // C_n^{n/2} => C20^10 (1)
    // nC2' (1)
    // i (1)
    // 2S_n^k for k=1..(n/2-1) => 2S20^k for k=1..9 (9)
    // σh (1)
    // nσv (1)
    // Total = 1+9+1+1+1+9+1+1 = 23.
    //
    // So there is NO extra separate "2σd" class here; only nσv (because for n>4 the common convention
    // doesn’t split v/d into distinct conjugacy classes). We'll use "20sv".

    std::vector<std::string> classNames;
    std::vector<int> classSizes;

    classNames.push_back("E");
    classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) {           // k=1..9
        classNames.push_back("2C20^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("C20^10");
    classSizes.push_back(1);

    classNames.push_back("20C2'");
    classSizes.push_back(n);

    classNames.push_back("i");
    classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) {           // k=1..9
        classNames.push_back("2S20^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("sh");
    classSizes.push_back(1);

    classNames.push_back("20sv");
    classSizes.push_back(n);

    auto make1D = [&](std::string name, int rSign, int sSign, int invSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(23);
 
        chi.push_back(1.0); // E
 
        for (int k = 1; k <= m - 1; ++k) {
            double val = (rSign == 1) ? 1.0 : ((k % 2) ? -1.0 : 1.0);
            chi.push_back(val);
        }
 
        double c10 = 1.0;                 // (-1)^10 = +1, so always +1 for r=-1 too
        chi.push_back(c10);               // C20^10
 
        chi.push_back(double(sSign));     // 20C2'
 
        chi.push_back(double(invSign));   // i
 
        for (int k = 1; k <= m - 1; ++k) {
            double val = (rSign == 1) ? 1.0 : ((k % 2) ? -1.0 : 1.0);
            chi.push_back(double(invSign) * val);   // 2S20^k
        }
 
        chi.push_back(double(invSign) * c10);       // sh
        chi.push_back(double(invSign) * double(sSign)); // 20sv
 
        return Ir{name, 1, chi};
    };

    auto make2D = [&](std::string name, int l, int invSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(23);
 
        auto rot = [&](int k) -> double {
            return 2.0 * std::cos(2.0 * pi * double(l) * double(k) / double(n));
        };
 
        chi.push_back(2.0); // E
 
        for (int k = 1; k <= m - 1; ++k) chi.push_back(rot(k));  // 2C20^k
 
        double c10 = 2.0 * ((l % 2) ? -1.0 : 1.0);               // 2(-1)^l
        chi.push_back(c10);                                      // C20^10
 
        chi.push_back(0.0);                                      // 20C2'
 
        chi.push_back(2.0 * double(invSign));                    // i
 
        for (int k = 1; k <= m - 1; ++k) chi.push_back(double(invSign) * rot(k)); // 2S20^k
 
        chi.push_back(double(invSign) * c10);                     // sh
        chi.push_back(0.0);                                       // 20sv
 
        return Ir{name, 2, chi};
    };


    std::vector<Ir> irreps;
    irreps.reserve(26);

    // 1D g
    irreps.push_back(make1D("A1g", +1, +1, +1));
    irreps.push_back(make1D("A2g", +1, -1, +1));
    irreps.push_back(make1D("B1g", -1, +1, +1));
    irreps.push_back(make1D("B2g", -1, -1, +1));

    // 2D g: E1g..E9g
    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l) + "g", l, +1));

    // 1D u
    irreps.push_back(make1D("A1u", +1, +1, -1));
    irreps.push_back(make1D("A2u", +1, -1, -1));
    irreps.push_back(make1D("B1u", -1, +1, -1));
    irreps.push_back(make1D("B2u", -1, -1, -1));

    // 2D u: E1u..E9u
    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l) + "u", l, -1));

    return PointGroupCharacterTable::build(
        "D20h",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_D40h()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // D40h: n=40, order = 4n = 160
    // Even n, standard class structure (non-splitting of C2' vs C2'' for n>4 in the convention
    // used for D20h above):
    //   E
    //   2C40^k, k=1..(n/2-1)=19
    //   C40^20
    //   40C2'
    //   i
    //   2S40^k, k=1..19
    //   sh
    //   40sv
    //
    // Total classes = 1 + 19 + 1 + 1 + 1 + 19 + 1 + 1 = 44
    constexpr int n = 40;
    constexpr int order = 4 * n;   // 160
    constexpr int m = n / 2;       // 20
    constexpr double pi = 3.14159265358979323846;

    std::vector<std::string> classNames;
    std::vector<int> classSizes;

    classNames.push_back("E");      classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) { // k=1..19
        classNames.push_back("2C40^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("C40^20"); classSizes.push_back(1);
    classNames.push_back("40C2'");  classSizes.push_back(n);

    classNames.push_back("i");      classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) { // k=1..19
        classNames.push_back("2S40^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("sh");     classSizes.push_back(1);
    classNames.push_back("40sv");   classSizes.push_back(n);

    // 1D reps determined by signs on generators:
    // r = C40, s = C2' (perpendicular 2-fold), inv = i.
    auto make1D = [&](std::string name, int rSign, int sSign, int invSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(44);

        auto rpow = [&](int k)->double {
            if (rSign == 1) return 1.0;
            return (k % 2) ? -1.0 : 1.0; // (-1)^k
        };

        chi.push_back(1.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rpow(k)); // 2C40^k

        double c20 = rpow(20);           // C40^20 (always +1 even if rSign=-1)
        chi.push_back(c20);

        chi.push_back(double(sSign));   // 40C2'

        chi.push_back(double(invSign)); // i

        for (int k = 1; k <= m - 1; ++k)
            chi.push_back(double(invSign) * rpow(k)); // 2S40^k = i*C40^k

        chi.push_back(double(invSign) * c20);          // sh = i*C40^20
        chi.push_back(double(invSign) * double(sSign)); // 40sv = i*C2'

        return Ir{name, 1, chi};
    };

    // 2D reps E_l, l=1..(m-1)=19:
    // χ(C40^k)=2cos(2π l k / 40), χ(C40^20)=2cos(π l)=2(-1)^l, χ(C2')=0
    // g/u: χ(i)=±2, χ(S40^k)=±χ(C40^k), χ(sh)=±χ(C40^20), χ(sv)=0
    auto make2D = [&](std::string name, int l, int invSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(44);

        auto rot = [&](int k)->double {
            return 2.0 * std::cos(2.0 * pi * double(l) * double(k) / double(n));
        };

        chi.push_back(2.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rot(k)); // 2C40^k

        double c20 = 2.0 * ((l % 2) ? -1.0 : 1.0); // 2(-1)^l
        chi.push_back(c20);                        // C40^20

        chi.push_back(0.0);                        // 40C2'

        chi.push_back(2.0 * double(invSign));      // i

        for (int k = 1; k <= m - 1; ++k)
            chi.push_back(double(invSign) * rot(k)); // 2S40^k

        chi.push_back(double(invSign) * c20);       // sh
        chi.push_back(0.0);                         // 40sv

        return Ir{name, 2, chi};
    };

    std::vector<Ir> irreps;
    irreps.reserve(8 + 2 * (m - 1)); // 8 one-dimensional + 38 two-dimensional = 46

    // 1D g
    irreps.push_back(make1D("A1g", +1, +1, +1));
    irreps.push_back(make1D("A2g", +1, -1, +1));
    irreps.push_back(make1D("B1g", -1, +1, +1));
    irreps.push_back(make1D("B2g", -1, -1, +1));

    // 2D g: E1g..E19g
    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l) + "g", l, +1));

    // 1D u
    irreps.push_back(make1D("A1u", +1, +1, -1));
    irreps.push_back(make1D("A2u", +1, -1, -1));
    irreps.push_back(make1D("B1u", -1, +1, -1));
    irreps.push_back(make1D("B2u", -1, -1, -1));

    // 2D u: E1u..E19u
    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l) + "u", l, -1));

    return PointGroupCharacterTable::build(
        "D40h",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_C3v()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "C3v", 6,
        {"E", "2C3", "3sv"},
        {  1,   2,     3  },
        { Ir{"A1", 1, { 1,  1,  1}},
          Ir{"A2", 1, { 1,  1, -1}},
          Ir{"E",  2, { 2, -1,  0}} }
    );
}
static PointGroupCharacterTable make_C4v()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // C4v: order = 8
    // Standard classes:
    //   E (1)
    //   2C4 (2)
    //   C2 (1)
    //   2σv (2)
    //   2σd (2)
    std::vector<std::string> classNames = {"E","2C4","C2","2sv","2sd"};
    std::vector<int> classSizes = {  1,   2,    1,   2,    2 };

    std::vector<Ir> irreps = {
        Ir{"A1",1,{ 1, 1, 1, 1, 1}},
        Ir{"A2",1,{ 1, 1, 1,-1,-1}},
        Ir{"B1",1,{ 1,-1, 1, 1,-1}},
        Ir{"B2",1,{ 1,-1, 1,-1, 1}},
        Ir{"E", 2,{ 2, 0,-2, 0, 0}},
    };

    return PointGroupCharacterTable::build(
        "C4v",
        8,
        classNames,
        classSizes,
        irreps
    );
}


static PointGroupCharacterTable make_C5v()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // C5v: order = 2n = 10
    // Classes:
    //   E (1)
    //   2C5 (2)        : C5, C5^4
    //   2C5^2 (2)      : C5^2, C5^3
    //   5σv (5)
    std::vector<std::string> classNames = {"E","2C5","2C5^2","5sv"};
    std::vector<int> classSizes        = {  1,   2,    2,    5 };

    // Irreps: A1, A2, E1, E2
    // For E_k (k=1,2): χ(C5^m) = 2cos(2π k m / 5), χ(σv)=0
    constexpr double c72  = 0.30901699437494745;  // cos(72°)
    constexpr double c144 = -0.8090169943749473;  // cos(144°)

    std::vector<Ir> irreps = {
        Ir{"A1",1,{ 1, 1, 1,  1}},
        Ir{"A2",1,{ 1, 1, 1, -1}},
        Ir{"E1",2,{ 2, 2*c72,  2*c144, 0}},   // {2, 0.6180339, -1.6180339, 0}
        Ir{"E2",2,{ 2, 2*c144, 2*c72,  0}},   // {2, -1.6180339, 0.6180339, 0}
    };

    return PointGroupCharacterTable::build(
        "C5v",
        10,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_C6v()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // C6v: order = 2n = 12
    // Standard classes (n even, not divisible by 4):
    //   E (1)
    //   2C6 (2)          (C6, C6^5)
    //   2C6^2 (2)        (C3, C3^2)
    //   C6^3 (1)         (C2)
    //   3σv (3)
    //   3σd (3)
    constexpr int n = 6;
    constexpr int order = 2 * n; // 12
    constexpr double pi = 3.14159265358979323846;

    std::vector<std::string> classNames = {
        "E","2C6","2C6^2","C6^3","3sv","3sd"
    };
    std::vector<int> classSizes = {
        1,  2,    2,     1,    3,    3
    };

    // 1D irreps: A1, A2, B1, B2
    // Let r=C6 with r->+1 for A, r->-1 for B; and mirror signs distinguish v/d:
    // In C6v the standard convention matches the usual table:
    //   A1: all +1
    //   A2: sv +1, sd -1
    //   B1: C2' like +1 pattern -> sv -1, sd +1 with C6->-1
    //   B2: sv -1, sd -1 with C6->-1 but differs on C6^2
    //
    // Easiest is to hardcode the standard C6v table for this class order:
    //  A1:  1  1  1  1  1  1
    //  A2:  1  1  1  1 -1 -1
    //  B1:  1 -1  1 -1  1 -1
    //  B2:  1 -1  1 -1 -1  1
    //  E1:  2  1 -1 -2  0  0
    //  E2:  2 -1 -1  2  0  0
    //
    // (E1 and E2 are the two 2D irreps.)

    std::vector<Ir> irreps = {
        Ir{"A1", 1, {1,  1,  1,  1,  1,  1}},
        Ir{"A2", 1, {1,  1,  1,  1, -1, -1}},
        Ir{"B1", 1, {1, -1,  1, -1,  1, -1}},
        Ir{"B2", 1, {1, -1,  1, -1, -1,  1}},
        Ir{"E1", 2, {2,  1, -1, -2,  0,  0}},
        Ir{"E2", 2, {2, -1, -1,  2,  0,  0}},
    };

    return PointGroupCharacterTable::build(
        "C6v",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_C8v()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // C8v: order = 2n = 16
    // For n divisible by 4, use split mirror classes:
    //   E
    //   2C8, 2C8^3
    //   C8^2, C8^4
    //   4σv, 4σd
    // Total classes = 1 + 1 + 1 + 1 + 1 + 1 + 1 = 7 classes? (sizes sum to 16)
    //
    // More explicitly (7 classes):
    //   E (1)
    //   2C8 (2)        [k=1]
    //   2C8^2 (2)      [k=2]
    //   2C8^3 (2)      [k=3]
    //   C8^4 (1)       [k=4]
    //   4sv (4)
    //   4sd (4)
    constexpr int n = 8;
    constexpr int order = 2 * n; // 16
    constexpr double pi = 3.14159265358979323846;

    std::vector<std::string> classNames = {
        "E","2C8","2C8^2","2C8^3","C8^4","4sv","4sd"
    };
    std::vector<int> classSizes = {
        1,  2,    2,     2,     1,   4,    4
    };

    // 1D irreps for Cnv with n divisible by 4: A1, A2, B1, B2
    // Use r = C8, and for mirrors distinguish σv vs σd.
    //
    // In 1D reps:
    //   r -> +1 (A) or -1 (B)
    //   σv -> +1 ("1") or -1 ("2")
    //   σd -> (r)^(n/4) * σv; for n=8, n/4=2 and (r)^2 is always +1 in 1D (even if r=-1),
    //         so σv and σd have the same character in these 1D reps.
    // (As with D8h, the σv/σd distinction shows up in 2D reps depending on convention/basis.)
    auto make1D = [&](std::string name, int rSign, int svSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(7);

        auto rpow = [&](int k)->double {
            if (rSign == 1) return 1.0;
            return (k % 2) ? -1.0 : 1.0; // (-1)^k
        };

        chi.push_back(1.0);      // E
        chi.push_back(rpow(1));  // 2C8
        chi.push_back(rpow(2));  // 2C8^2
        chi.push_back(rpow(3));  // 2C8^3
        chi.push_back(rpow(4));  // C8^4

        chi.push_back(double(svSign)); // 4sv
        chi.push_back(double(svSign)); // 4sd (same in these 1D reps)

        return Ir{name, 1, chi};
    };

    // 2D irreps: E_l, l=1..(n/2-1)=3
    // χ(C8^k)=2cos(2π l k / 8), χ(C8^4)=2(-1)^l
    // For reflections in Cnv 2D irreps: χ(σ) = 0.
    auto make2D = [&](std::string name, int l) -> Ir {
        std::vector<double> chi;
        chi.reserve(7);

        auto rot = [&](int k)->double {
            return 2.0 * std::cos(2.0 * pi * double(l) * double(k) / double(n));
        };

        chi.push_back(2.0);     // E
        chi.push_back(rot(1));  // 2C8
        chi.push_back(rot(2));  // 2C8^2
        chi.push_back(rot(3));  // 2C8^3
        double c4 = 2.0 * ((l % 2) ? -1.0 : 1.0); // 2(-1)^l
        chi.push_back(c4);      // C8^4

        chi.push_back(0.0);     // 4sv
        chi.push_back(0.0);     // 4sd

        return Ir{name, 2, chi};
    };

    std::vector<Ir> irreps;
    irreps.reserve(4 + 3); // 7 irreps total

    irreps.push_back(make1D("A1", +1, +1));
    irreps.push_back(make1D("A2", +1, -1));
    irreps.push_back(make1D("B1", -1, +1));
    irreps.push_back(make1D("B2", -1, -1));

    for (int l = 1; l <= 3; ++l)
        irreps.push_back(make2D("E" + std::to_string(l), l));

    return PointGroupCharacterTable::build(
        "C8v",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_C10v()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // C10v: order = 2n = 20
    // Even n class structure:
    //   E
    //   2C10^k, k=1..(n/2-1)=4
    //   C10^5
    //   10σv
    // Total classes = 1 + 4 + 1 + 1 = 7
    constexpr int n = 10;
    constexpr int order = 2 * n;     // 20
    constexpr int m = n / 2;         // 5
    constexpr double pi = 3.14159265358979323846;

    std::vector<std::string> classNames;
    std::vector<int> classSizes;

    classNames.push_back("E");
    classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) { // k=1..4
        classNames.push_back("2C10^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("C10^5");
    classSizes.push_back(1);

    classNames.push_back("10sv");
    classSizes.push_back(n);

    // 1D irreps (n even): A1, A2, B1, B2
    // Take r = C10 and s = σv. Then:
    //   r -> +1 (A) or -1 (B)
    //   s -> +1 ("1") or -1 ("2")
    auto make1D = [&](std::string name, int rSign, int sSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(7);

        auto rpow = [&](int k)->double {
            if (rSign == 1) return 1.0;
            return (k % 2) ? -1.0 : 1.0; // (-1)^k
        };

        chi.push_back(1.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rpow(k)); // 2C10^k

        chi.push_back(rpow(5));          // C10^5 ( = -1 if r=-1 and k odd; here k=5 odd)
        chi.push_back(double(sSign));    // 10sv

        return Ir{name, 1, chi};
    };

    // 2D irreps: E_l, l=1..(m-1)=4
    // χ(C10^k)=2cos(2π l k / 10), χ(C10^5)=2(-1)^l, χ(σv)=0
    auto make2D = [&](std::string name, int l) -> Ir {
        std::vector<double> chi;
        chi.reserve(7);

        auto rot = [&](int k)->double {
            return 2.0 * std::cos(2.0 * pi * double(l) * double(k) / double(n));
        };

        chi.push_back(2.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rot(k)); // 2C10^k

        double c5 = 2.0 * ((l % 2) ? -1.0 : 1.0); // 2(-1)^l
        chi.push_back(c5);

        chi.push_back(0.0); // 10sv

        return Ir{name, 2, chi};
    };

    std::vector<Ir> irreps;
    irreps.reserve(4 + (m - 1)); // 4 one-dimensional + 4 two-dimensional = 8 irreps

    irreps.push_back(make1D("A1", +1, +1));
    irreps.push_back(make1D("A2", +1, -1));
    irreps.push_back(make1D("B1", -1, +1));
    irreps.push_back(make1D("B2", -1, -1));

    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l), l));

    return PointGroupCharacterTable::build(
        "C10v",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_C20v()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // C20v: order = 2n = 40
    // Even n class structure:
    //   E
    //   2C20^k, k=1..(n/2-1)=9
    //   C20^10
    //   20σv
    // Total classes = 1 + 9 + 1 + 1 = 12
    constexpr int n = 20;
    constexpr int order = 2 * n;     // 40
    constexpr int m = n / 2;         // 10
    constexpr double pi = 3.14159265358979323846;

    std::vector<std::string> classNames;
    std::vector<int> classSizes;

    classNames.push_back("E");
    classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) { // k=1..9
        classNames.push_back("2C20^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("C20^10");
    classSizes.push_back(1);

    classNames.push_back("20sv");
    classSizes.push_back(n);

    // 1D irreps (n even): A1, A2, B1, B2
    // Take r = C20 and s = σv. Then for 1D reps:
    //   r -> +1 (A) or -1 (B)
    //   s -> +1 ("1") or -1 ("2")
    auto make1D = [&](std::string name, int rSign, int sSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(12);

        auto rpow = [&](int k)->double {
            if (rSign == 1) return 1.0;
            return (k % 2) ? -1.0 : 1.0; // (-1)^k
        };

        chi.push_back(1.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rpow(k)); // 2C20^k

        chi.push_back(rpow(10));         // C20^10 (always +1 even if r=-1)
        chi.push_back(double(sSign));    // 20sv

        return Ir{name, 1, chi};
    };

    // 2D irreps: E_l, l=1..(m-1)=9
    // χ(C20^k) = 2cos(2π l k / 20), χ(C20^10)=2(-1)^l, χ(σv)=0
    auto make2D = [&](std::string name, int l) -> Ir {
        std::vector<double> chi;
        chi.reserve(12);

        auto rot = [&](int k)->double {
            return 2.0 * std::cos(2.0 * pi * double(l) * double(k) / double(n));
        };

        chi.push_back(2.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rot(k)); // 2C20^k

        double c10 = 2.0 * ((l % 2) ? -1.0 : 1.0); // 2(-1)^l
        chi.push_back(c10);

        chi.push_back(0.0); // 20sv

        return Ir{name, 2, chi};
    };

    std::vector<Ir> irreps;
    irreps.reserve(4 + (m - 1)); // 4 one-dimensional + 9 two-dimensional = 13 irreps

    irreps.push_back(make1D("A1", +1, +1));
    irreps.push_back(make1D("A2", +1, -1));
    irreps.push_back(make1D("B1", -1, +1));
    irreps.push_back(make1D("B2", -1, -1));

    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l), l));

    return PointGroupCharacterTable::build(
        "C20v",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_C40v()
{
    using Ir = PointGroupCharacterTable::Irrep;

    // C40v: order = 2n = 80
    // For even n, standard class structure:
    //   E
    //   2C40^k, k=1..(n/2-1)=19
    //   C40^20
    //   40σv
    // Total classes = 1 + 19 + 1 + 1 = 22
    constexpr int n = 40;
    constexpr int order = 2 * n;     // 80
    constexpr int m = n / 2;         // 20
    constexpr double pi = 3.14159265358979323846;

    std::vector<std::string> classNames;
    std::vector<int> classSizes;

    classNames.push_back("E");
    classSizes.push_back(1);

    for (int k = 1; k <= m - 1; ++k) { // k=1..19
        classNames.push_back("2C40^" + std::to_string(k));
        classSizes.push_back(2);
    }

    classNames.push_back("C40^20");
    classSizes.push_back(1);

    classNames.push_back("40sv");
    classSizes.push_back(n);

    // 1D reps: A1, A2, B1, B2 (since n even)
    // Use r = C40, s = σv. For 1D reps:
    //   r -> +1 (A) or -1 (B)
    //   s -> +1 ("1") or -1 ("2")
    // Then χ(C40^k)=(r)^k, χ(σv)=s
    auto make1D = [&](std::string name, int rSign, int sSign) -> Ir {
        std::vector<double> chi;
        chi.reserve(22);

        auto rpow = [&](int k)->double {
            if (rSign == 1) return 1.0;
            return (k % 2) ? -1.0 : 1.0; // (-1)^k
        };

        chi.push_back(1.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rpow(k)); // 2C40^k

        chi.push_back(rpow(20));          // C40^20 (always +1 even if r=-1)
        chi.push_back(double(sSign));     // 40sv

        return Ir{name, 1, chi};
    };

    // 2D reps: E_l, l=1..(m-1)=19
    // χ(C40^k)=2cos(2π l k / 40), χ(C40^20)=2(-1)^l, χ(σv)=0
    auto make2D = [&](std::string name, int l) -> Ir {
        std::vector<double> chi;
        chi.reserve(22);

        auto rot = [&](int k)->double {
            return 2.0 * std::cos(2.0 * pi * double(l) * double(k) / double(n));
        };

        chi.push_back(2.0); // E

        for (int k = 1; k <= m - 1; ++k) chi.push_back(rot(k)); // 2C40^k

        double c20 = 2.0 * ((l % 2) ? -1.0 : 1.0); // 2(-1)^l
        chi.push_back(c20);

        chi.push_back(0.0); // 40sv

        return Ir{name, 2, chi};
    };

    std::vector<Ir> irreps;
    irreps.reserve(4 + (m - 1)); // 4 one-dimensional + 19 two-dimensional = 23 irreps

    irreps.push_back(make1D("A1", +1, +1));
    irreps.push_back(make1D("A2", +1, -1));
    irreps.push_back(make1D("B1", -1, +1));
    irreps.push_back(make1D("B2", -1, -1));

    for (int l = 1; l <= m - 1; ++l)
        irreps.push_back(make2D("E" + std::to_string(l), l));

    return PointGroupCharacterTable::build(
        "C40v",
        order,
        classNames,
        classSizes,
        irreps
    );
}

static PointGroupCharacterTable make_D3h()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "D3h", 12,
        {"E","2C3","3C2","sh","2S3","3sv"},
        {  1,  2,   3,   1,   2,   3 },
        { Ir{"A1'",  1, { 1,  1,  1,  1,  1,  1}},
          Ir{"A2'",  1, { 1,  1, -1,  1,  1, -1}},
          Ir{"E'",   2, { 2, -1,  0,  2, -1,  0}},
          Ir{"A1''", 1, { 1,  1,  1, -1, -1, -1}},
          Ir{"A2''", 1, { 1,  1, -1, -1, -1,  1}},
          Ir{"E''",  2, { 2, -1,  0, -2,  1,  0}} }
    );
}

static PointGroupCharacterTable make_Td()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "Td", 24,
        // Classes: E, 8C3, 3C2, 6S4, 6σd
        {"E", "8C3", "3C2", "6S4", "6sd"},
        {  1,    8,    3,    6,    6 },
        { Ir{"A1", 1, {1,  1,  1,  1,  1}},
          Ir{"A2", 1, {1,  1,  1, -1, -1}},
          Ir{"E",  2, {2, -1,  2,  0,  0}},
          Ir{"T1", 3, {3,  0, -1,  1, -1}},
          Ir{"T2", 3, {3,  0, -1, -1,  1}} }
    );
}

static PointGroupCharacterTable make_Oh()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "Oh", 48,
        // Classes: E, 8C3, 6C2, 6C4, 3C2(=C4^2), i, 6S4, 8S6, 3sh, 6sd
        {"E","8C3","6C2","6C4","3C2","i","6S4","8S6","3sh","6sd"},
        {  1,   8,    6,    6,    3,   1,   6,    8,    3,    6 },
        { Ir{"A1g",1,{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
          Ir{"A2g",1,{ 1, 1,-1,-1, 1, 1,-1, 1, 1,-1}},
          Ir{"Eg", 2,{ 2,-1, 0, 0, 2, 2, 0,-1, 2, 0}},
          Ir{"T1g",3,{ 3, 0,-1, 1,-1, 3, 1, 0,-1,-1}},
          Ir{"T2g",3,{ 3, 0, 1,-1,-1, 3,-1, 0,-1, 1}},
          Ir{"A1u",1,{ 1, 1, 1, 1, 1,-1,-1,-1,-1,-1}},
          Ir{"A2u",1,{ 1, 1,-1,-1, 1,-1, 1,-1,-1, 1}},
          Ir{"Eu", 2,{ 2,-1, 0, 0, 2,-2, 0, 1,-2, 0}},
          Ir{"T1u",3,{ 3, 0,-1, 1,-1,-3,-1, 0, 1, 1}},
          Ir{"T2u",3,{ 3, 0, 1,-1,-1,-3, 1, 0, 1,-1}} }
    );
}

static PointGroupCharacterTable make_Th()
{
    using Ir = PointGroupCharacterTable::Irrep;
    return PointGroupCharacterTable::build(
        "Th", 24,
        {"E","8C3","3C2","i","8S6","3sh"},
        { 1,   8,    3,   1,   8,    3 },
        {
          Ir{"Ag",1,{ 1, 1, 1, 1, 1, 1}},
          Ir{"Au",1,{ 1, 1, 1,-1,-1,-1}},
          Ir{"Eg",2,{ 2,-1, 2, 2,-1, 2}},
          Ir{"Eu",2,{ 2,-1, 2,-2, 1,-2}},
          Ir{"Tg",3,{ 3, 0,-1, 3, 0,-1}},
          Ir{"Tu",3,{ 3, 0,-1,-3, 0, 1}}
        }
    );
}

static PointGroupCharacterTable make_Ih()
{
    using Ir = PointGroupCharacterTable::Irrep;

    const double phi  = (1.0 + std::sqrt(5.0)) / 2.0;
    const double phi2 = (1.0 - std::sqrt(5.0)) / 2.0;

    return PointGroupCharacterTable::build(
        "Ih", 120,
        {"E","12C5","12C5^2","20C3","15C2","i","12S10","12S10^3","20S6","15sh"},
        { 1,   12,     12,     20,    15,   1,   12,      12,      20,    15 },

        {
          Ir{"Ag",1,{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
          Ir{"Au",1,{ 1, 1, 1, 1, 1,-1,-1,-1,-1,-1}},

          Ir{"T1g",3,{ 3, phi, phi2, 0,-1, 3, phi, phi2, 0,-1}},
          Ir{"T1u",3,{ 3, phi, phi2, 0,-1,-3,-phi,-phi2, 0, 1}},

          Ir{"T2g",3,{ 3, phi2, phi, 0,-1, 3, phi2, phi, 0,-1}},
          Ir{"T2u",3,{ 3, phi2, phi, 0,-1,-3,-phi2,-phi, 0, 1}},

          Ir{"Gg",4,{ 4,-1,-1, 1, 0, 4,-1,-1, 1, 0}},
          Ir{"Gu",4,{ 4,-1,-1, 1, 0,-4, 1, 1,-1, 0}},

          Ir{"Hg",5,{ 5, 0, 0,-1, 1, 5, 0, 0,-1, 1}},
          Ir{"Hu",5,{ 5, 0, 0,-1, 1,-5, 0, 0, 1,-1}}
        }
    );
}

PointGroupCharacterTable PointGroupCharacterTable::from_symbol(const std::string& symbol)
{
    if (symbol == "C1")  return make_C1();
    if (symbol == "Cs")  return make_Cs();
    if (symbol == "Ci")  return make_Ci();
    if (symbol == "C2")  return make_C2();
    if (symbol == "C2v") return make_C2v();
    if (symbol == "D2h") return make_D2h();
    if (symbol == "D4h") return make_D4h();
    if (symbol == "D6h") return make_D6h();
    if (symbol == "D8h") return make_D8h();
    if (symbol == "D10h") return make_D10h();
    if (symbol == "D20h") return make_D20h();
    if (symbol == "D40h") return make_D40h();
    if (symbol == "C3v") return make_C3v();
    if (symbol == "C4v") return make_C4v();
    if (symbol == "C5v") return make_C5v();
    if (symbol == "C6v") return make_C6v();
    if (symbol == "C8v") return make_C8v();
    if (symbol == "C10v") return make_C10v();
    if (symbol == "C20v") return make_C20v();
    if (symbol == "C40v") return make_C40v();
    if (symbol == "D3h") return make_D3h();
    if (symbol == "Td")  return make_Td();
    if (symbol == "Oh")  return make_Oh();
    if (symbol == "Th")  return make_Th();
    if (symbol == "Ih")  return make_Ih();

    std::cerr << " Warning character table not tabulated: unsupported point group '" << symbol 
          << "' — falling back to C1\n";
    return make_C1();

    throw std::runtime_error("PointGroupCharacterTable: unsupported group: " + symbol);
}

} // namespace pwdft
