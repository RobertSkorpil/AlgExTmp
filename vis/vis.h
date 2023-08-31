#pragma once
#include <cstdlib>
#include <array>

namespace vis
{
    using vec = std::array<double, 4>;

    struct point
    {
        vec pos;
        vec color;
    };

    void init();
    void draw_point(const point& p);
    void draw_points(const point* p, size_t n);
    void swap();
}
