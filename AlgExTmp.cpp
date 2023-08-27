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
            return -1.0;
        else
            return 1.0;
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

template<symb::Expr M>
constexpr auto inverse2(M m)
{
    using namespace symb;
    auto p2{ perm_sign2() };
    auto det_g{ Sum<'a', 2>(Sum<'b', 2>(Ix<'a'>(Ix<'b'>(p2)) * Ix<'a'>(I<0>(m)) * Ix<'b'>(I<1>(m)))) };
    auto i{ Inv(det_g) };
    return simplify(arr(arr(i * I<1>(I<1>(m)), i * -I<1>(I<0>(m))), arr(i * -I<0>(I<1>(m)), i * I<0>(I<0>(m)))));
}

int main()
{
    SetConsoleOutputCP(CP_UTF8);

    using namespace symb;

    auto vars{ arr(v<'r'>, v<'a'>) };
    std::cout << to_str(vars.to_str()) << '\n';

    auto p2{ perm_sign2() };
    std::cout << to_str(p2.to_str()) << '\n';

    auto t{ For<'x', 2>(For<'b', 2>(Sum<'a', 2>(d(Ix<'a'>(Ix<'b'>(p2))) / d(Ix<'x'>(vars))))) };
    std::cout << to_str(t.to_str()) << '\n';

    auto p3{ perm_sign3() };
    std::cout << to_str(p3.to_str()) << '\n';

    //    auto vars2{ arr(arr(v<'a'>, v<'b'>, v<'c'>), arr(v<'d'>, v<'e'>, v<'f'>), arr(v<'g'>, v<'h'>, v<'i'>)) };

    auto x{ d(Inv(v<'x'>)) / d(v<'x'>) };
    std::cout << to_str(x.to_str()) << '\n';

    constexpr auto D{ vars.size() };

    auto map{ arr(v<'r'> *Cos(v<'a'>), v<'r'> *Sin(v<'a'>)) };
    std::cout << to_str(map.to_str()) << '\n';

    auto J{ For<'v', D>(For<'u', D>(d(Ix<'v'>(map)) / d(Ix<'u'>(vars)))) };
    std::cout << to_str(J.to_str()) << '\n';

    auto g{ For<'a', D>(For<'b', D>(Sum<'c', D>((Ix<'a'>(Ix<'c'>(J)) * Ix<'b'>(Ix<'c'>(J)))))) };
    std::cout << to_str(g.to_str()) << '\n';

    auto ginv{ inverse2(g) };
    std::cout << to_str(ginv.to_str()) << '\n';


    double R{}, A{};
    struct eval_context
    {
        double &r;
        double &a;
        std::array<double, 4> out;

        eval_context(double& r, double& a)
            : r{ r }, a{ a }
        {
        }

        auto var_value(const var& v, std::optional<size_t> index) const
        {
            switch (v.n)
            {
            case 'r':
                return r;
            case 'a':
                return a;
            }
        }

        void set_var_value(const var& v, std::optional<size_t> index, scalar value) const
        {
        }
    };
    eval_context ctx{ R, A };

    auto Gamma{ For<'i', D>(For<'k', D>(For<'l', D>(Sum<'m', D>(
        c<0.5> *Ix<'i'>(Ix<'m'>(ginv)) *
    (d(Ix<'k'>(Ix<'m'>(g))) / d(Ix<'l'>(vars)) + d(Ix<'l'>(Ix<'m'>(g))) / d(Ix<'k'>(vars)) - d(Ix<'k'>(Ix<'l'>(g))) / d(Ix<'m'>(vars))))))) };
    std::cout << to_str(Gamma.to_str()) << '\n';


    double vr{}, va{ 0.1 };
    double ar, aa;



    std::array<std::array<std::array<double, D>, D>, D> gamma_val;
    gamma_val[0][0][0] = I<0>(I<0>(I<0>(Gamma))).eval(ctx);
    gamma_val[0][0][1] = I<0>(I<0>(I<1>(Gamma))).eval(ctx);
    gamma_val[0][1][0] = I<0>(I<1>(I<0>(Gamma))).eval(ctx);
    gamma_val[0][1][1] = I<0>(I<1>(I<1>(Gamma))).eval(ctx);
    gamma_val[1][0][0] = I<1>(I<0>(I<0>(Gamma))).eval(ctx);
    gamma_val[1][0][1] = I<1>(I<0>(I<1>(Gamma))).eval(ctx);
    gamma_val[1][1][0] = I<1>(I<1>(I<0>(Gamma))).eval(ctx);
    gamma_val[1][1][1] = I<1>(I<1>(I<1>(Gamma))).eval(ctx);
}
