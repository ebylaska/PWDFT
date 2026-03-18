#include "PointGroupCharacterTables.hpp"
#include <utility> // std::move

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

PointGroupCharacterTable PointGroupCharacterTable::from_symbol(const std::string& symbol)
{
    if (symbol == "C1")  return make_C1();
    if (symbol == "Cs")  return make_Cs();
    if (symbol == "Ci")  return make_Ci();
    if (symbol == "C2")  return make_C2();
    if (symbol == "C2v") return make_C2v();
    if (symbol == "D2h") return make_D2h();
    if (symbol == "C3v") return make_C3v();
    if (symbol == "D3h") return make_D3h();
    if (symbol == "Td")  return make_Td();
    if (symbol == "Oh")  return make_Oh();

    throw std::runtime_error("PointGroupCharacterTable: unsupported group: " + symbol);
}

} // namespace pwdft
