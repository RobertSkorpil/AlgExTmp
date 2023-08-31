#pragma once
#include "symb.h"
#include <iostream>

template<size_t... I>
constexpr double perm_sign()
{
    constexpr std::array<size_t, sizeof...(I)> ix {{ I... }};

    double sign{ 1 };
    for (size_t i{}; i < ix.size(); ++i)
        for (size_t j{ i + 1 }; j < ix.size(); ++j)
            if (ix[j] > ix[i])
                sign *= -1;
            else if (ix[j] == ix[i])
                return 0;

    return sign;
}

template<size_t D, size_t... I>
constexpr auto perm_tensor();

template<size_t D, size_t i, size_t... I>
constexpr auto perm_tensor_arr()
{
    using namespace symb;
    if constexpr (i == D)
        return nothing{};
    else
        return Array<decltype(perm_tensor<D, I..., i>()), decltype(perm_tensor_arr<D, i + 1, I...>())>{};
}

template<size_t D, size_t... I>
constexpr auto perm_tensor()
{
    using namespace symb;
    if constexpr (sizeof...(I) == D)
        return Constant<perm_sign<I...>()>{};
    else
        return perm_tensor_arr<D, 0, I...>();
}

template<size_t i, size_t D, symb::Expr PermT>
constexpr auto det_indexed_perm(PermT perm)
{
    if constexpr (i == D)
        return perm;
    else
        return Ix<'a' + i>(det_indexed_perm<i + 1, D>(perm));
}

template<size_t i, symb::Expr MatrixT, symb::Expr PermT>
constexpr auto det_indexed_factor(MatrixT matrix, PermT perm);

template<size_t i, symb::Expr MatrixT, symb::Expr PermT>
constexpr auto det_indexed_factor(MatrixT matrix, PermT perm)
{
    constexpr auto D{ MatrixT::size() };
    using namespace symb;

    if constexpr (i == D)
        return det_indexed_perm<0, D>(perm);
    else
        return Ix<'a' + i>(I<i>(matrix)) * det_indexed_factor<i + 1>(matrix, perm);
}

template<size_t i, symb::Expr MatrixT, symb::Expr FactorT>
constexpr auto det_sum(MatrixT matrix, FactorT factor);

template<size_t i, symb::Expr MatrixT, symb::Expr FactorT>
constexpr auto det_sum(MatrixT matrix, FactorT factor)
{
    constexpr auto D{ MatrixT::size() };
    using namespace symb;

    if constexpr (i == D)
        return factor;
    else
        return Sum<'a' + i, D>(det_sum<i + 1>(matrix, factor));
}

template<symb::Expr MatrixT>
constexpr auto matrix_determinant(MatrixT matrix)
{
    using namespace symb;
    constexpr auto D{ MatrixT::size() };

    if constexpr (D == 1)
        return I<0>(I<0>(matrix));
    else
    {
        auto perm{ perm_tensor<MatrixT::size()>() };
        auto factor{ det_indexed_factor<0>(matrix, perm) };
        return det_sum<0>(matrix, factor);
    }
}

template<size_t i, size_t j, symb::Expr MatrixT>
constexpr auto mminor(MatrixT matrix)
{
    using namespace symb;
    constexpr auto D{ MatrixT::size() };
    return ForSkip<'i', D, i>(ForSkip<'j', D, j>(Ix<'j'>(Ix<'i'>(matrix))));
}

template<size_t i, size_t j, size_t D, symb::Expr MatrixT, symb::Expr DetT>
constexpr auto matrix_inverse2(MatrixT matrix, DetT det_i);

template<size_t i, size_t j, size_t D, symb::Expr MatrixT, symb::Expr DetT>
constexpr auto matrix_inverse2(MatrixT matrix, DetT det_i)
{
    using namespace symb;
    if constexpr (j == D)
        return nothing{};
    else
    {
        constexpr auto f{ (i + j) % 2 ? -1.0 : 1.0 };
        return Array<decltype(det_i * c<f> * matrix_determinant(mminor<j, i>(matrix))), decltype(matrix_inverse2<i, j + 1, D>(matrix, det_i))>{};
    }
}

template<size_t i, size_t D, symb::Expr MatrixT, symb::Expr DetT>
constexpr auto matrix_inverse1(MatrixT matrix, DetT det_i);

template<size_t i, size_t D, symb::Expr MatrixT, symb::Expr DetT>
constexpr auto matrix_inverse1(MatrixT matrix, DetT det_i)
{
    using namespace symb;
    if constexpr (i == D)
        return nothing{};
    else
        return Array<decltype(matrix_inverse2<i, 0, D>(matrix, det_i)), decltype(matrix_inverse1<i + 1, D>(matrix, det_i))>{};
}

template<symb::Expr MatrixT>
constexpr auto matrix_inverse(MatrixT matrix)
{
    using namespace symb;
    auto det_i{ simplify(Inv(matrix_determinant(matrix))) };

    constexpr auto D{ MatrixT::size() };
    return simplify(matrix_inverse1<0, D>(matrix, det_i));
}
