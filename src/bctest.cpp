#include "tools/bcconv.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>

double testfun(const double xi, const double yi, double& fx, double& fy)
{
  // evaluate function x**3 + y**3 at xi, yi
  fx = 3.*std::pow(xi,2.);
  fy = 3.*std::pow(yi,2.);
  return std::pow(xi,3.) + std::pow(yi,3.);
}

int main(int argc, char** argv)
{

  // set grid size
  const int nx = 100;
  const int ny = 100;

  // set grid resolution
  const double gridres = 0.01;

  // set up data grid
  std::vector<std::vector<double> > grid(ny,std::vector<double>(nx,0.));
  std::vector<std::vector<double> > gridx(ny,std::vector<double>(nx,0.));
  std::vector<std::vector<double> > gridy(ny,std::vector<double>(nx,0.));
  for (int iy=0; iy<ny; ++iy)
    for (int ix=0; ix<nx; ++ix)
    {
      double sx, sy;
      grid[iy][ix]  = testfun(ix*gridres, iy*gridres, sx, sy);
      gridx[iy][ix] = sx;
      gridy[iy][ix] = sy;
    }

  // evaluate cubic spline at some grid points
  for (int iy=0; iy<ny; iy+=10)
  {
    for (int ix=0; ix<nx; ix+=10)
    {
      // get x,y values
      const double x = static_cast<double>(ix)*gridres;
      const double y = static_cast<double>(iy)*gridres;
      // get true function value
      const double fval  = grid[iy][ix];
      const double fxval = gridx[iy][ix];
      const double fyval = gridy[iy][ix];
      // evaluate cubic spline at x,y
      double sval, sx, sy, sxx, syy, sxy;
      if (!bicubicConvolution(x, y, grid, gridres, sval, sx, sy, sxx, syy, sxy))
        continue;

      // std::cout.precision(4);
      // std::cout.width(10);
      std::cout << "x,y: "
        << std::setprecision(3) << std::setw(5) << x << ", "
        << std::setprecision(3) << std::setw(5) << y << " = "
        << " fval: " << std::setprecision(3) << std::setw(5) << fval
        << " sval: " << std::setprecision(3) << std::setw(5) << sval
        << " fx/sx: " << fyval << " / " << sxy
        << std::endl;

    }
  }

    // performance test
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> dist(1, nx-3);
    const int N = 10000;
    auto start = std::chrono::system_clock::now();
    for (int i=0; i<N; ++i)
    {
      const double x = static_cast<double>(dist(mt))*gridres;
      const double y = static_cast<double>(dist(mt))*gridres;
      double sval, sx, sy, sxx, syy, sxy;
      bicubicConvolution(x, y, grid, gridres, sval, sx, sy, sxx, syy, sxy);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << "Average call took: " << diff.count()/static_cast<double>(N)
              << " seconds" << std::endl;
}
