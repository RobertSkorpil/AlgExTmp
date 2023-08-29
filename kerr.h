#pragma once
#include <array>
#include <functional>

using vec = std::array<double, 4>;
using point = std::array<double, 4>;
using tensor3 = std::array<std::array<std::array<double, 4>, 4>, 4>;
struct kerr
{
  std::function<double(vec v, point p)> magnitude;
  std::function<tensor3(point p)> christoffel;

  kerr(double M, double J);
};
