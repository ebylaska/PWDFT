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


static double det3(const SymOp& op)
{
    const auto& R = op.R;
    return R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])
         - R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])
         + R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);
}

static void mul3(const double A[3][3], const double B[3][3], double C[3][3])
{
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++) {
            C[i][j]=0.0;
            for(int k=0;k<3;k++) C[i][j]+=A[i][k]*B[k][j];
        }
}

static bool is_I3(const double R[3][3], double eps)
{
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++) {
            double target = (i==j)?1.0:0.0;
            if (std::fabs(R[i][j]-target) > eps) return false;
        }
    return true;
}

// Td class index order must match your table: {"E","8C3","3C2","6S4","6sd"}
static int td_class_index(const SymOp& op)
{
    const double eps = 1e-6;
    const double (&R)[3][3] = op.R;

    if (is_I3(R, eps)) return 0;

    const double d  = det3(op);
    const double tr = R[0][0] + R[1][1] + R[2][2];

    if (d > 0.0) {
        if (std::fabs(tr - 0.0) < 1e-4) return 1; // C3
        if (std::fabs(tr + 1.0) < 1e-4) return 2; // C2
    } else {
        double R2[3][3];
        mul3(R, R, R2);
        if (is_I3(R2, 1e-5)) return 4; // reflection σd
        return 3;                      // S4
    }

    return -1;
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

    // Build class labels consistent with the character table
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


/*******************************************
 *                                         *
 *         Symmetry::rotate_operators      *
 *                                         *
 *******************************************/
void Symmetry::rotate_operators(const double *U)
{
    double UT[9];

    // transpose U
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            UT[3*i+j] = U[3*j+i];

    for(auto &op : ops_)
    {
        double R0[9], tmp[9], R1[9];

        // copy operator matrix
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                R0[3*i+j] = op.R[i][j];

        // tmp = U * R0
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            {
                tmp[3*i+j] = 0.0;
                for(int k=0;k<3;k++)
                    tmp[3*i+j] += U[3*i+k] * R0[3*k+j];
            }

        // R1 = tmp * U^T
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            {
                R1[3*i+j] = 0.0;
                for(int k=0;k<3;k++)
                    R1[3*i+j] += tmp[3*i+k] * UT[3*k+j];
            }

        // overwrite operator
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                op.R[i][j] = R1[3*i+j];
    }
}



} //namespace pwdft
