#pragma once

#include "symb.h"
#include "matrix.h"
#include "affine_space.h"

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

template<size_t D>
struct riemann
{
    using vec = affine_space<D>::vec;
    using point = affine_space<D>::point;
    using tensor2 = affine_space<D>::tensor2;
    using tensor3 = affine_space<D>::tensor3;

    static vec affine_acceleration(tensor3 gamma, vec v)
    {
        vec a{};
        for (int i{}; i < D; ++i)
            for (int j{}; j < D; ++j)
                for (int k{}; k < D; ++k)
                    a[i] += (-gamma[i][j][k]) * v[j] * v[k];
        return a;
    }

    static vec coordinate_acceleration(tensor3 gamma, vec v)
    {
        vec a{};
        for (int i{}; i < D; ++i)
            for (int j{}; j < D; ++j)
                for (int k{}; k < D; ++k)
                    a[i] += (-gamma[i][j][k] + gamma[0][j][k] * v[i]) * v[j] * v[k];
        return a;
    }
};
