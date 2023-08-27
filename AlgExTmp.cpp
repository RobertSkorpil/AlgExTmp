#include <iostream>
#include <locale>
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#include <array>
#include "symb.h"
#include "str.h"

template<size_t i, size_t j, size_t k>
constexpr auto perm_sign3_3()
{
    using namespace symb;
    constexpr double sign{ [&]() constexpr {
        double s = -1.0;
        if constexpr (i == j || j == k || i == k)
            return 0.0;
        if constexpr (j < i)
            s *= -1;
        if constexpr (k < j)
            s *= -1;
        if constexpr (k < i)
            s *= -1;
        return s;
    }() };

    return Constant<sign>{};
}

template<size_t i, size_t j>
constexpr auto perm_sign3_2()
{
    using namespace symb;
    return arr(perm_sign3_3<i, j, 0>(), perm_sign3_3<i, j, 1>(), perm_sign3_3<i, j, 2>());
}

template<size_t i>
constexpr auto perm_sign3_1()
{
    using namespace symb;
    return arr(perm_sign3_2<i, 0>(), perm_sign3_2<i, 1>(), perm_sign3_2<i, 2>());
}

constexpr auto perm_sign3()
{
    using namespace symb;
    return arr(perm_sign3_1<0>(), perm_sign3_1<1>(), perm_sign3_1<2>());
}

template<size_t i, size_t j>
constexpr auto perm_sign2_2()
{
    using namespace symb;
    constexpr double sign{ [&]() constexpr {
        if constexpr (i == j)
            return 0.0;
        if constexpr (j < i)
            return 1.0;
        else
            return -1.0;
    }() };

    return Constant<sign>{};
}

template<size_t i>
constexpr auto perm_sign2_1()
{
    using namespace symb;
    return arr(perm_sign2_2<i, 0>(), perm_sign2_2<i, 1>());
}

constexpr auto perm_sign2()
{
    using namespace symb;
    return arr(perm_sign2_1<0>(), perm_sign2_1<1>());
}

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
constexpr auto cofactor(MatrixT matrix)
{
    using namespace symb;
    constexpr auto D{ MatrixT::size() };

    return ForSkip<'i', D, i>(ForSkip<'j', D, j>(Ix<'i'>(Ix<'j'>(matrix))));
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
        return Array<decltype(det_i *  matrix_determinant(cofactor<i, j>(matrix))), decltype(matrix_inverse2<i, j + 1, D>(matrix, det_i))>{};
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

template<symb::Expr M>
constexpr auto inverse2(M m)
{
    using namespace symb;
    auto p2{ perm_sign2() };
    auto det_g{ Sum<'a', 2>(Sum<'b', 2>(Ix<'a'>(Ix<'b'>(p2)) * Ix<'a'>(I<0>(m)) * Ix<'b'>(I<1>(m)))) };
    auto i{ Inv(det_g) };
    return simplify(arr(arr(i * I<1>(I<1>(m)), i * -I<1>(I<0>(m))), arr(i * -I<0>(I<1>(m)), i * I<0>(I<0>(m)))));
}

template<symb::Expr MetricT, symb::Expr Vars>
constexpr auto Christoffel(MetricT g, Vars vars)
{
    using namespace symb;
    constexpr auto D{ g.size() };
    auto ginv{ inverse2(g) };
    return For<'i', D>(For<'k', D>(For<'l', D>(Sum<'m', D>(
        c<0.5> * Ix<'i'>(Ix<'m'>(ginv)) *
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

int main()
{
    SetConsoleOutputCP(CP_UTF8);

    using namespace symb;
    
//    auto matrix{ arr(arr(v<'a'>, v<'b'>, v<'c'>, v<'d'>, v<'e'>), arr(v<'f'>, v<'g'>, v<'h'>, v<'i'>, v<'j'>), arr(v<'k'>, v<'l'>, v<'m'>, v<'n'>, v<'o'>), arr(v<'p'>, v<'q'>, v<'r'>, v<'s'>, v<'t'>), arr(v<'u'>, v<'v'>, v<'w'>, v<'x'>, v<'y'>)) };
    auto matrix{ arr(arr(v<'a'>, v<'b'>, v<'c'>, v<'d'>), arr(v<'e'>, v<'f'>, v<'g'>, v<'h'>), arr(v<'i'>, v<'j'>, v<'k'>, v<'l'>), arr(v<'m'>, v<'n'>, v<'o'>, v<'p'>)) };
//    auto matrix{ arr(arr(v<'a'>, c<0.0>, c<0.0>, c<0.0>), arr(c<0.0>, v<'f'>, c<0.0>, c<0.0>), arr(c<0.0>, c<0.0>, v<'k'>, c<0.0>), arr(c<0.0>, c<0.0>, c<0.0>, v<'p'>)) };
//    auto matrix{ arr(arr(v<'a'>, v<'b'>, v<'c'>), arr(v<'d'>, v<'e'>, v<'f'>), arr(v<'g'>, v<'h'>, v<'i'>)) };
//    auto matrix{ arr(arr(v<'a'>, v<'b'>), arr(v<'c'>, v<'d'>)) };
    std::cout << to_str(matrix.to_str()) << '\n';

//    std::cout << to_str(perm_tensor<4>().to_str()) << '\n';
//    std::cout << to_str(cofactor<2, 2>(matrix).to_str()) << '\n';
    
    auto f{ matrix_inverse(matrix) };
    std::cout << to_str(f.to_str()) << '\n';
#if 0
#if 0
    auto t{ arr(
        arr(arr(v<'a'>, v<'b'>, v<'c'>), arr(v<'d'>, v<'e'>, v<'f'>), arr(v<'g'>, v<'h'>, v<'i'>)),
        arr(arr(v<'j'>, v<'k'>, v<'l'>), arr(v<'m'>, v<'n'>, v<'o'>), arr(v<'p'>, v<'q'>, v<'r'>)),
        arr(arr(v<'s'>, v<'t'>, v<'u'>), arr(v<'v'>, v<'w'>, v<'x'>), arr(v<'y'>, v<'z'>, v<'Z'>))
    )
        };
    auto t2{ For<'a', 3>(For<'b', 3>(For<'c', 3>(For<'d', 3>(Sum<'l', 3>(Ix<'b'>(Ix<'d'>(Ix<'l'>(t)))))))) };
    auto R{ Riemann(t, c<0.0>) };

    auto tst{ simplify(Sin(v<'x'>) * Cos(v<'x'>)) };
    std::cout << to_str(tst.to_str()) << '\n';
#endif

    auto vars{ arr(v<'r'>, v<'a'>) };
    std::cout << to_str(vars.to_str()) << '\n';

    auto p2{ perm_sign2() };
    std::cout << to_str(p2.to_str()) << '\n';

    auto t{ For<'x', 2>(For<'b', 2>(Sum<'a', 2>(d(Ix<'a'>(Ix<'b'>(p2))) / d(Ix<'x'>(vars))))) };
    std::cout << to_str(t.to_str()) << '\n';

    auto p3{ perm_sign3() };
    std::cout << to_str(p3.to_str()) << '\n';

    auto vars2{ arr(arr(v<'a'>, v<'b'>, v<'c'>), arr(v<'d'>, v<'e'>, v<'f'>), arr(v<'g'>, v<'h'>, v<'i'>)) };
    auto cof1{ ForSkip<'a', 3, 1>(ForSkip<'b', 3, 1>(Ix<'b'>(Ix<'a'>(vars2)))) };
    std::cout << to_str(cof1.to_str()) << '\n';

    auto x{ d(Inv(v<'x'>)) / d(v<'x'>) };
    std::cout << to_str(x.to_str()) << '\n';

    constexpr auto D{ vars.size() };

    auto map{ arr(v<'r'> *Cos(v<'a'>), v<'r'> *Sin(v<'a'>)) };
    std::cout << to_str(map.to_str()) << '\n';

    auto J{ For<'v', D>(For<'u', D>(d(Ix<'v'>(map)) / d(Ix<'u'>(vars)))) };
    std::cout << to_str(J.to_str()) << '\n';

    auto g{ For<'a', D>(For<'b', D>(Sum<'c', D>((Ix<'a'>(Ix<'c'>(J)) * Ix<'b'>(Ix<'c'>(J)))))) };
    std::cout << to_str(g.to_str()) << '\n';

    auto vars_sphere{ arr(v<'a'>, v<'b'>) };
    auto g_sphere{ arr(arr(v<'R'> *v<'R'>, c<0.0>), arr(c<0.0>, v<'R'> *v<'R'> *Sin(v<'a'>) * Sin(v<'a'>))) };
    auto gamma_sphere{ Christoffel(g_sphere, vars_sphere) };
    std::cout << to_str(gamma_sphere.to_str()) << '\n';

    auto riemann_sphere{ Riemann(gamma_sphere, vars_sphere) };
    std::cout << to_str(riemann_sphere.to_str()) << '\n';

    double R{ 1.0 }, A{ 0.5 }, B{};
    struct eval_context
    {
        double& r;
        double& a;
        double& b;
        std::array<double, 4> out{};

        eval_context(double& r, double& a, double& b)
            : r{ r }, a{ a }, b{ b }
        {
        }

        auto var_value(const var& v, std::optional<size_t> index) const
        {
            switch (v.n)
            {
            case 'R':
                return r;
            case 'a':
                return a;
            }
        }

        void set_var_value(const var& v, std::optional<size_t> index, scalar value) const
        {
        }
    };
    eval_context ctx{ R, A, B };

    auto Gamma{ Christoffel(g, vars) };
    std::cout << to_str(Gamma.to_str()) << '\n';

    auto riemann_polar{ Riemann(Gamma, vars) };
    std::cout << to_str(riemann_polar.to_str()) << '\n';


    double v[2]{ 0.0, 0.2 };
    double a[2];

    auto G{ gamma_sphere };

    for (int p{}; p < 10000000; p++)
    {
        double tdelta = 1e-4;
        std::array<std::array<std::array<double, D>, D>, D> gamma_val;
        gamma_val[0][0][0] = I<0>(I<0>(I<0>(G))).eval(ctx);
        gamma_val[1][0][0] = I<0>(I<0>(I<1>(G))).eval(ctx);
        gamma_val[0][1][0] = I<0>(I<1>(I<0>(G))).eval(ctx);
        gamma_val[1][1][0] = I<0>(I<1>(I<1>(G))).eval(ctx);
        gamma_val[0][0][1] = I<1>(I<0>(I<0>(G))).eval(ctx);
        gamma_val[1][0][1] = I<1>(I<0>(I<1>(G))).eval(ctx);
        gamma_val[0][1][1] = I<1>(I<1>(I<0>(G))).eval(ctx);
        gamma_val[1][1][1] = I<1>(I<1>(I<1>(G))).eval(ctx);

        a[0] = 0;
        a[1] = 0;
        for (int i{}; i < D; ++i)
            for (int j{}; j < D; ++j)
                for (int k{}; k < D; ++k)
                    a[i] += -gamma_val[i][j][k] * v[j] * v[k];

        A += tdelta * v[0];
        B += tdelta * v[1];

        v[0] += tdelta * a[0];
        v[1] += tdelta * a[1];

        auto x{ R * cos(A) };
        auto y{ R * sin(A) };
        if(p % 100 == 0)
            std::cout << "A: " << A << " B: " << B << '\n';
    }
#endif
}
