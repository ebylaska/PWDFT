#ifndef PHASE_CLASSIFY_MINIMAL_HPP
#define PHASE_CLASSIFY_MINIMAL_HPP

#include <cmath>
#include <queue>
#include <vector>

namespace pwdft {

inline double phase_det3(const double a[9])
{
    return a[0]*(a[4]*a[8] - a[5]*a[7])
         - a[1]*(a[3]*a[8] - a[5]*a[6])
         + a[2]*(a[3]*a[7] - a[4]*a[6]);
}

inline bool phase_invert3x3(const double a[9], double inv[9])
{
    double d = phase_det3(a);
    if (std::abs(d) < 1e-14) return false;

    inv[0] =  (a[4]*a[8] - a[5]*a[7]) / d;
    inv[1] = -(a[1]*a[8] - a[2]*a[7]) / d;
    inv[2] =  (a[1]*a[5] - a[2]*a[4]) / d;

    inv[3] = -(a[3]*a[8] - a[5]*a[6]) / d;
    inv[4] =  (a[0]*a[8] - a[2]*a[6]) / d;
    inv[5] = -(a[0]*a[5] - a[2]*a[3]) / d;

    inv[6] =  (a[3]*a[7] - a[4]*a[6]) / d;
    inv[7] = -(a[0]*a[7] - a[1]*a[6]) / d;
    inv[8] =  (a[0]*a[4] - a[1]*a[3]) / d;

    return true;
}

inline double phase_pbc_dist(const double* r, int i, int j, const double cell[9], const double inv[9])
{
    double dx = r[3*i]   - r[3*j];
    double dy = r[3*i+1] - r[3*j+1];
    double dz = r[3*i+2] - r[3*j+2];

    double sx = inv[0]*dx + inv[1]*dy + inv[2]*dz;
    double sy = inv[3]*dx + inv[4]*dy + inv[5]*dz;
    double sz = inv[6]*dx + inv[7]*dy + inv[8]*dz;

    sx -= std::round(sx);
    sy -= std::round(sy);
    sz -= std::round(sz);

    dx = cell[0]*sx + cell[3]*sy + cell[6]*sz;
    dy = cell[1]*sx + cell[4]*sy + cell[7]*sz;
    dz = cell[2]*sx + cell[5]*sy + cell[8]*sz;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

inline const char* classify_minimal(int nion, const double* rion, const double* unita)
{
    if (nion <= 0 || rion == nullptr || unita == nullptr) {
        return "gas";
    }

    double cell[9] = {
        unita[0], unita[1], unita[2],
        unita[3], unita[4], unita[5],
        unita[6], unita[7], unita[8]
    };

    double omega = std::abs(phase_det3(cell));
    double rho = (omega > 1e-12) ? static_cast<double>(nion) / omega : 0.0;

    if (nion == 1) {
        return (rho < 0.01) ? "molecule_in_box" : "condensed";
    }

    double inv[9];
    if (!phase_invert3x3(cell, inv)) {
        return "condensed";
    }

    std::vector<double> dmin(nion, 1e100);
    for (int i = 0; i < nion; ++i) {
        for (int j = i + 1; j < nion; ++j) {
            double d = phase_pbc_dist(rion, i, j, cell, inv);
            if (d < dmin[i]) dmin[i] = d;
            if (d < dmin[j]) dmin[j] = d;
        }
    }

    double mean_dmin = 0.0;
    for (double d : dmin) mean_dmin += d;
    mean_dmin /= static_cast<double>(nion);

    double bond_cut = 1.3 * mean_dmin;

    std::vector<std::vector<int>> g(nion);
    for (int i = 0; i < nion; ++i) {
        for (int j = i + 1; j < nion; ++j) {
            double d = phase_pbc_dist(rion, i, j, cell, inv);
            if (d < bond_cut) {
                g[i].push_back(j);
                g[j].push_back(i);
            }
        }
    }

    std::vector<int> vis(nion, 0);
    int ncomp = 0;
    int largest = 0;

    for (int i = 0; i < nion; ++i) {
        if (vis[i]) continue;

        ++ncomp;
        int sz = 0;
        std::queue<int> q;
        q.push(i);
        vis[i] = 1;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            ++sz;

            for (int v : g[u]) {
                if (!vis[v]) {
                    vis[v] = 1;
                    q.push(v);
                }
            }
        }

        if (sz > largest) largest = sz;
    }

    if (rho < 0.01) {
        if (ncomp == 1 && largest == nion) return "molecule_in_box";
        return "gas";
    }

    return "condensed";
}

} // namespace pwdft

#endif
