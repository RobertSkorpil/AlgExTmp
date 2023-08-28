#include "kerr.h"
#ifndef __INTELLISENSE__
#include "symb.h"
#endif
#include "riemann.h"
#include <iostream>

using namespace symb;
namespace
{
    struct context
    {
        double J, M, t, r, theta, phi;

        static constexpr auto _2{ c<2.0> };

        double var_value(const var& v, std::optional<size_t> index) const
        {
            switch (v.n)
            {
            case 'J':
                return J;
            case 'M':
                return M;
            case 't':
                return t;
            case 'r':
                return r;
            case u'Θ':
                return theta;
            case u'φ':
                return phi;
            default:
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
    };
}

std::function<std::array<std::array<std::array<double, 4>, 4>, 4>(double t, double r, double theta, double phi)>
kerr(double M_, double J_)
{
    auto t{ v<'t'> };
    auto r{ v<'r'> };
    auto theta{ v<u'Θ'> };
    auto phi{ v<u'φ'> };

    auto J{ v<'J'> };
    auto M{ v<'M'> };
    auto _0{ c<0.0> };
    auto _1{ c<1.0> };
    auto _2{ c<2.0> };
    auto a{ J / (_2 * M) };
    auto S{ r * r + a * a * Cos(theta) * Cos(theta) };
    auto D{ r * r - _2 * M * r + a * a };
    auto g{ arr(
        arr(-(_1 - _2 * M * r / S), _0, _0, -_2 * M * r * a * Sin(theta) * Sin(theta) / S),
        arr(_0, S / D, _0, _0),
        arr(_0, _0, S, _0),
        arr(-_2 * M * r * a * Sin(theta) * Sin(theta) / S, _0, _0, (r * r + a * a + _2 * M * r * a * a * Sin(theta) * Sin(theta) / S) * Sin(theta) * Sin(theta))
    ) };

    auto vars{ arr(t, r, theta, phi) };
    auto G{ Christoffel(g, vars) };

    std::cout << "Metric: \n" << to_str(g.to_str()) << '\n';
    std::cout << "Γ:      \n" << to_str(G.to_str()) << '\n';

    return [=](double t_, double r_, double theta_, double phi_)
    {
        context ctx{ J_, M_, t_, r_, theta_, phi_ };
        std::array<std::array<std::array<double, 4>, 4>, 4> val;

        G.eval_array(ctx, &val[0][0][0]);

        return val;
    };
}
