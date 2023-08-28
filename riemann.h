#pragma once

#include "symb.h"
#include "matrix.h"

template<symb::Expr MetricT, symb::Expr Vars>
constexpr auto Christoffel(MetricT g, Vars vars)
{
    using namespace symb;
    constexpr auto D{ g.size() };
    auto ginv{ matrix_inverse(g) };
    return For<'i', D>(For<'k', D>(For<'l', D>(Sum<'m', D>(
        c<1.0> / c<2.0> * Ix<'i'>(Ix<'m'>(ginv)) *
    (d(Ix<'k'>(Ix<'m'>(g))) / d(Ix<'l'>(vars)) + d(Ix<'l'>(Ix<'m'>(g))) / d(Ix<'k'>(vars)) - d(Ix<'l'>(Ix<'k'>(g))) / d(Ix<'m'>(vars)))))));
}

template<symb::Expr ChristoffelT, symb::Expr Vars>
constexpr auto Riemann(ChristoffelT G, Vars vars)
{
    using namespace symb;
    constexpr auto D{ G.size() };

    return For<'a', D>(For<'b', D>(For<'c', D>(For<'d', D>(
        d(Ix<'b'>(Ix<'d'>(Ix<'a'>(G)))) / d(Ix<'c'>(vars))
      - d(Ix<'b'>(Ix<'c'>(Ix<'a'>(G)))) / d(Ix<'d'>(vars))
      + Sum<'l', D>(Ix<'l'>(Ix<'c'>(Ix<'a'>(G))) * Ix<'b'>(Ix<'d'>(Ix<'l'>(G))))
      - Sum<'l', D>(Ix<'l'>(Ix<'d'>(Ix<'a'>(G))) * Ix<'b'>(Ix<'c'>(Ix<'l'>(G))))
    ))));
}

