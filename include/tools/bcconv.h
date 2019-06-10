///< Bicubic convolutional algorithm
///< https://en.wikipedia.org/wiki/Bicubic_interpolation

#pragma once

#include <cmath>
#include <vector>

/**
 * @brief evaluates cubic convolutional algorithm on 1d data
 * @param t spline coordinate (between 0 and 1)
 * @param f function values at grid points (-1,0,1,2)
 * @param s resulting spline value at t
 * @param sp first derivative of spline
 * @param spp second derivative of spline
 */
inline void cubicConvolution(const double t,
                             const double* f,
                             double& s,
                             double& sp,
                             double& spp)
{
  // get necessary powers of t
  const double t2 = std::pow(t,2);
  const double t3 = std::pow(t,3);

  const double k0 = +2.*f[1];
  const double k1 = -1.*f[0] + 1.*f[2];
  const double k2 = +2.*f[0] - 5.*f[1] + 4.*f[2] - 1.*f[3];
  const double k3 = -1.*f[0] + 3.*f[1] - 3.*f[2] + 1.*f[3];

  // get spline value at t and first two derivatives
  s   = 0.5 * (   k0 +    k1*t +    k2*t2 + k3*t3 );
  sp  = 0.5 * (   k1 + 2.*k2*t + 3.*k3*t2);
  spp = 0.5 * (2.*k2 + 6.*k3*t);
}

/**
 * @brief evaluates bicubic convolutional algorithm at given point x,y
 * @param x x-coordinate of evaluation point
 * @param y y-coordinate of evaluation point
 * @param grid regular grid containing data
 * @param resolution resolution of grid in x/y direction
 * @param s spline value at x,y
 */
inline bool bicubicConvolution(const double x,
                               const double y,
                               const std::vector<std::vector<double> >& grid,
                               const double resolution,
                               double& s,
                               double& sx,
                               double& sy,
                               double& sxx,
                               double& syy,
                               double& sxy)
{
  // compute indices into grid for given x,y coordinates
  const int xindex = static_cast<int>(x / resolution);
  const int yindex = static_cast<int>(y / resolution);

  // compute spline values tx,ty within grid
  const double tx = (x - static_cast<double>(xindex)*resolution) / resolution;
  const double ty = (y - static_cast<double>(yindex)*resolution) / resolution;

  // make sure that both indices valid
  if (xindex < 1 || yindex < 1)
    return false;
  if (yindex > static_cast<int>(grid.size())-3 ||
      xindex > static_cast<int>(grid[0].size())-3)
    return false;

  // evaluate 1d cubic convolutional algorithm
  std::vector<double> b(4,0.);
  std::vector<double> bp(4,0.);
  std::vector<double> bpp(4,0.);
  cubicConvolution(tx, &grid[yindex-1].data()[xindex-1], b[0], bp[0], bpp[0]);
  cubicConvolution(tx, &grid[yindex+0].data()[xindex-1], b[1], bp[1], bpp[1]);
  cubicConvolution(tx, &grid[yindex+1].data()[xindex-1], b[2], bp[2], bpp[2]);
  cubicConvolution(tx, &grid[yindex+2].data()[xindex-1], b[3], bp[3], bpp[3]);

  // evaluation in y-direction
  double sxxy, sxyy, sxxyy;
  cubicConvolution(ty, b.data(),   s,   sy,   syy);
  cubicConvolution(ty, bp.data(),  sx,  sxy,  sxyy);
  cubicConvolution(ty, bpp.data(), sxx, sxxy, sxxyy);

  sx /= resolution;
  sy /= resolution;
  sxx /= std::pow(resolution,2.);
  syy /= std::pow(resolution,2.);
  sxy /= std::pow(resolution,2.);

  return true;

}
