#ifndef UTIL_ASCII_PLOT_HPP
#define UTIL_ASCII_PLOT_HPP

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace pwdft {
namespace util_ascii {

/******************************************
 *                                        *
 *          _util_ascii_plotter_          *
 *                                        *
 ******************************************/
/*
 * Draw a simple ASCII plot of y(x) to the selected output stream.
 *
 * The plotting window is controlled by util_ascii_setwindow().
 * Data are mapped into a fixed MaxRow x MaxCol character grid.
 *
 * Notes:
 *   - The first character of 'symbol' is used as the plotting mark.
 *   - The plot window is shared through namespace-local static state,
 *     replacing the original Fortran common block.
 */

inline constexpr int MaxRow = 20;
inline constexpr int MaxCol = 86;

/* plotting window, replacing the old Fortran common block */
inline double XL = 0.0;
inline double XR = 1.0;
inline double YD = 0.0;
inline double YU = 1.0;

/******************************************
 *                                        *
 *         _util_ascii_setwindow_         *
 *                                        *
 ******************************************/
/*
 * Set the plotting window limits used by the ASCII plot routines.
 */
inline void util_ascii_setwindow(double xmin, double xmax,
                                 double ymin, double ymax)
{
    XL = xmin;
    XR = xmax;
    YD = ymin;
    YU = ymax;
}

/******************************************
 *                                        *
 *          _util_ascii_xscale_           *
 *                                        *
 ******************************************/
/*
 * Map an x value from plot coordinates to the ASCII plot column index.
 * The returned index follows the original Fortran 1-based convention.
 */
inline int ascii_xscale(double x)
{
    return static_cast<int>(
        std::llround((x - XL) * (MaxCol - 1) / (XR - XL) + 1.0)
    );
}

/******************************************
 *                                        *
 *          _util_ascii_yscale_           *
 *                                        *
 ******************************************/
/*
 * Map a y value from plot coordinates to the ASCII plot row index.
 * The returned index follows the original Fortran 1-based convention.
 */
inline int ascii_yscale(double y)
{
    return static_cast<int>(
        std::llround((y - YU) * (MaxRow - 1) / (YD - YU) + 1.0)
    );
}

/* helper matching Fortran E10.3-like formatting */
inline std::string format_e10_3(double x)
{
    std::ostringstream os;
    os << std::uppercase << std::scientific << std::setprecision(3) << x;
    return os.str();
}

/******************************************
 *                                        *
 *          _util_ascii_plotter_          *
 *                                        *
 ******************************************/
/*
 * Plot n points (x,y) as an ASCII graph.
 *
 * Inputs:
 *   mark   - prefix written at the start of each output line
 *   os     - output stream
 *   n      - number of points
 *   x, y   - input data arrays
 *   symbol - plotting symbol; only symbol[0] is used
 *   title  - plot title
 *   xlabel - x-axis label
 *   ylabel - y-axis label
 */
inline void util_ascii_plotter(const std::string& mark,
                               std::ostream& os,
                               int n,
                               const double* x,
                               const double* y,
                               const std::string& symbol,
                               const std::string& title,
                               const std::string& xlabel,
                               const std::string& ylabel)
{
    std::vector<std::string> point(MaxRow, std::string(MaxCol, ' '));

    const char psym = symbol.empty() ? '*' : symbol[0];

    auto center_pad = [](const std::string& s, int width) -> int {
        int l = static_cast<int>(s.size());
        int pad = width - (l / 2);
        if (pad < 0) pad = 0;
        return pad;
    };

    const int ltitle = center_pad(title, 43);
    const int lxlab  = center_pad(xlabel, 43);

    const std::string nstru = format_e10_3(YU);
    const std::string nstrd = format_e10_3(YD);

    /* set y-axis */
    {
        int col = std::min(MaxCol, 11 + ascii_xscale(0.0));
        col = std::max(1, col);
        for (int j = 1; j <= MaxRow; ++j)
            point[j - 1][col - 1] = ':';
    }

    /* set x-axis */
    for (int i = 12; i <= MaxCol; ++i) {
        int j = std::min(MaxRow, ascii_yscale(0.0));
        if (j >= 1)
            point[j - 1][i - 1] = '-';
    }

    /* set ylabels */
    for (int i = 1; i <= 10; ++i) {
        if (i <= static_cast<int>(nstru.size())) point[0][i - 1] = nstru[i - 1];
        if (i <= static_cast<int>(nstrd.size())) point[MaxRow - 1][i - 1] = nstrd[i - 1];
    }

    if (ascii_yscale(0.0) < 20) {
        int i = std::min(MaxRow, ascii_yscale(0.0));
        if (i >= 1) {
            point[i - 1][0] = ' ';
            point[i - 1][1] = '0';
            point[i - 1][2] = '.';
            point[i - 1][3] = '0';
            point[i - 1][4] = '0';
            point[i - 1][5] = '0';
            point[i - 1][6] = 'E';
            point[i - 1][7] = '+';
            point[i - 1][8] = '0';
            point[i - 1][9] = '0';
        }
    }

    /* plot points */
    for (int i = 0; i < n; ++i) {
        int row = std::min(ascii_yscale(y[i]), MaxRow);
        int col = std::min(11 + ascii_xscale(x[i]), MaxCol);
        row = std::max(1, row);
        col = std::max(1, col);
        point[row - 1][col - 1] = psym;
    }

    /* write graph */
    os << mark << std::string(ltitle, ' ') << title << '\n';
    os << mark << ylabel << '\n';

    for (int i = 0; i < MaxRow; ++i)
        os << mark << point[i] << '\n';

    os << mark << "           "
       << "|....................................|"
       << "....................................|" << '\n';

    {
        std::ostringstream line;
        line << mark
             << std::string(5, ' ')
             << std::uppercase << std::scientific << std::setprecision(3)
             << std::setw(10) << XL
             << std::string(27, ' ')
             << std::setw(10) << (XL + XR) / 2.0
             << std::string(27, ' ')
             << std::setw(10) << XR;
        os << line.str() << '\n';
    }

    os << mark << std::string(lxlab, ' ') << xlabel << '\n';
}

/******************************************
 *                                        *
 *      _util_ascii_plotter_vector_       *
 *                                        *
 ******************************************/
/*
 * Convenience overload using std::vector<double>.
 */
inline void util_ascii_plotter(const std::string& mark,
                               std::ostream& os,
                               const std::vector<double>& x,
                               const std::vector<double>& y,
                               const std::string& symbol,
                               const std::string& title,
                               const std::string& xlabel,
                               const std::string& ylabel)
{
    const int n = static_cast<int>(std::min(x.size(), y.size()));
    util_ascii_plotter(mark, os, n, x.data(), y.data(),
                       symbol, title, xlabel, ylabel);
}

} // namespace util_ascii
} // namespace pwdft

#endif
