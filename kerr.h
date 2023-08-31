#pragma once
#include <array>
#include <functional>
#include "affine_space.h"

struct kerr
{
    using vec = affine_space<4>::vec;
    using point = affine_space<4>::point;
    using tensor3 = affine_space<4>::tensor3;

    std::function<double(vec v, point p)> magnitude;
    std::function<tensor3(point p)> christoffel;
    std::function<point(point)> sphere_to_cart;
    std::function<double(vec)> ds;

    kerr(double M, double J);
};
