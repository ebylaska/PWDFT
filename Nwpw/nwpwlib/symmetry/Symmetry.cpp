#include <fstream>
#include <sstream>
#include <stdexcept>
#include <array>
#include <cmath>


#include "SpaceGroupDB.hpp"
#include "Symmetry.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

static bool has_identity(const std::vector<SymOp>& ops)
{
    const double tol = 1e-12;
    for (const auto& op : ops)
    {
        bool ok = true;
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                if (std::abs(op.R[i][j] - (i==j)) > tol)
                    ok = false;
        for (int i=0;i<3;i++)
            if (std::abs(op.t[i]) > tol)
                ok = false;
        if (ok) return true;
    }
    return false;
}


/*******************************************
 *                                         *
 *           Symmetry::Symmetry            *
 *                                         *
 *******************************************/
Symmetry::Symmetry()
{
    name_ = "identity";
    type_ = Type::Trivial;
    ita_number_ = 0;
    num_centering_ = 1;

    ops_.resize(1);
    auto& op = ops_[0];

    // Identity rotation
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            op.R[i][j] = (i == j);

    // Zero translation
    op.t[0] = op.t[1] = op.t[2] = 0.0;
}



/*******************************************
 *                                         *
 *           Symmetry::Symmetry            *
 *                                         *
 *******************************************/
Symmetry::Symmetry(const std::string& sg_name)
{
    const auto& db = spacegroup_db();
    const auto& matches = db.by_name(sg_name);

    if (matches.empty())
        throw std::runtime_error("Unknown space group: " + sg_name);

    const SpaceGroup& sg = *matches.front();

    type_ = Type::SpaceGroup;
    name_ = sg.name;
    ita_number_ = sg.number;
    num_centering_ = sg.num_centering;
    ops_ = sg.ops;

    if (!has_identity(ops_))
    {
       SymOp id;
       for (int i=0;i<3;i++)
       {
          for (int j=0;j<3;j++)
             id.R[i][j] = (i==j);
          id.t[i] = 0.0;
       }
       ops_.insert(ops_.begin(), id);
    }
}


/*******************************************
 *                                         *
 *           Symmetry::Symmetry            *
 *                                         *
 *******************************************/
 Symmetry::Symmetry(int ita_number)
{
    const auto& db = spacegroup_db();
    const auto& matches = db.by_number(ita_number);

    if (matches.empty())
        throw std::runtime_error(
            "Unknown space group number: " + std::to_string(ita_number)
        );

    const SpaceGroup& sg = *matches.front();

    type_ = Type::SpaceGroup;
    name_ = sg.name;
    ita_number_ = sg.number;
    num_centering_ = sg.num_centering;
    ops_ = sg.ops;

    if (!has_identity(ops_))
    {
       SymOp id;
       for (int i=0;i<3;i++)
       {
          for (int j=0;j<3;j++)
             id.R[i][j] = (i==j);
          id.t[i] = 0.0;
       }
       ops_.insert(ops_.begin(), id);
    }
}





/*******************************************
 *                                         *
 *        Symmetry::from_point_group       *
 *                                         *
 *******************************************/
Symmetry Symmetry::from_point_group(const std::string& pg_symbol)
{
    Symmetry s;
    s.type_ = Type::PointGroup;
    s.name_ = pg_symbol;
    s.ops_ = PointGroupGenerators::generate(pg_symbol);
    s.num_centering_ = 1;
    return s;
}

/*******************************************
 *                                         *
 *        Symmetry::from_json              *
 *                                         *
 *******************************************/
Symmetry Symmetry::from_json(const nlohmann::json& j)
{
    Symmetry s;

    // Source tells us intent
    const std::string source = j.value("source", "identity");

    if (source == "identity")
    {
        return Symmetry();  // Trivial
    }

    // If group info exists, prefer database construction
    if (j.contains("group_number"))
    {
        return Symmetry(j["group_number"].get<int>());
    }

    if (j.contains("group_name"))
    {
        return Symmetry(j["group_name"].get<std::string>());
    }

    // Otherwise build directly from operators
    s.type_ = Type::SpaceGroup;
    s.name_ = j.value("name", "custom");
    s.ops_.clear();

    for (const auto& jop : j.at("ops"))
    {
        SymOp op;
        for (int i = 0; i < 3; ++i)
           for (int j = 0; j < 3; ++j)
              op.R[i][j] = jop["R"][i][j].get<double>();

        op.t[0] = jop["t"][0].get<double>();
        op.t[1] = jop["t"][1].get<double>();
        op.t[2] = jop["t"][2].get<double>();

        s.ops_.push_back(op);
    }

    s.num_centering_ = 1; // recompute later if needed
    return s;
}




/*******************************************
 *                                         *
 *            Symmetry::order              *
 *                                         *
 *******************************************/
size_t Symmetry::order() const
{
    return ops_.size();
}

/*******************************************
 *                                         *
 *            Symmetry::operators          *
 *                                         *
 *******************************************/
const std::vector<SymOp>& Symmetry::operators() const
{
    return ops_;
}

/*******************************************
 *                                         *
 *            Symmetry::name               *
 *                                         *
 *******************************************/
const std::string& Symmetry::name() const
{
    return name_;
}

/*******************************************
 *                                         *
 *       Symmetry::is_space_group          *
 *                                         *
 *******************************************/
bool Symmetry::is_space_group() const
{
    return type_ == Type::SpaceGroup;
}

/*******************************************
 *                                         *
 *       Symmetry::is_point_group          *
 *                                         *
 *******************************************/
bool Symmetry::is_point_group() const
{
    return type_ == Type::PointGroup;
}

/*******************************************
 *                                         *
 *         Symmetry::num_centering         *
 *                                         *
 *******************************************/
int Symmetry::num_centering() const
{
    return num_centering_;
}



} //namespace pwdft
