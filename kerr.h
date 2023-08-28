#pragma once
#include <array>
#include <functional>

std::function<std::array<std::array<std::array<double, 4>, 4>, 4>(double t, double r, double theta, double phi)>
kerr(double M, double J);
