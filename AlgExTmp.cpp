#include <iostream>
#include <locale>
#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#endif
#include <array>
#include "symb.h"
#include "matrix.h"
#include "riemann.h"
#include "kerr.h"
#include "str.h"

//#define MATRIX_INVERSE
//#define MATRIX_INVERSE_4
//#define SPHERE_SURFACE
//#define SCHWARZSCHILD
//#define SCHWARZSCHILD_RIEMANN
//#define POLAR
#define KERR


namespace {
#ifdef SCHWARZSCHILD
    void schwarzschild()
    {
        using namespace symb;
        auto one{ c<1.0> };
        auto nul{ c<0.0> };
        auto C{ v<'c'> };
        auto t{ v<'t'> };
        auto r{ v<'r'> };
        auto theta{ v<u'Θ'> };
        auto phi{ v<u'φ'> };

        auto vars{ arr(t, r, theta, phi) };

        auto rs{ c<2.0> *v<'M'> };

        auto g{ arr(arr(-(one - rs / r), nul, nul, nul), arr(nul, one / (one - rs / r), nul, nul), arr(nul, nul, r * r, nul), arr(nul, nul, nul, r * r * Sin(theta) * Sin(theta))) };
        auto G{ Christoffel(g, vars) };

        std::cout << "Metric: \n" << to_str(g.to_str()) << '\n';
        std::cout << "Γ:      \n" << to_str(G.to_str()) << '\n';

#ifdef SCHWARZSCHILD_RIEMANN
        auto R{ Riemann(G, vars) };
        std::cout << "R:      \n" << to_str(R.to_str()) << '\n';
#endif
    }
#endif

#ifdef KERR
    void kerr_geodesic()
    {
        double J{ 1.999 };
        double M{ 1. };
        kerr k { M, J };

        std::array<double, 4> x { {0., 2.7, 3.1415926536 / 2, 0 } };
        std::array<double, 4> v { {1.0, 0., 0., 0.1615} }; // coordinate velocity (v[0] = dt/dt = 1)

        size_t cnt{};
        for (;;)
        {
            double tdelta{ 1e-4 };
            double a[4]{ {} };

            auto GV{ k.christoffel(x) };
            for (int i{}; i < 4; ++i)
                for (int j{}; j < 4; ++j)
                    for (int k{}; k < 4; ++k)
                        a[i] += (-GV[i][j][k] + GV[0][j][k] * v[i]) * v[j] * v[k];

            for (size_t i{}; i < 4; ++i)
                x[i] += tdelta * v[i];

            for (size_t i{}; i < 4; ++i)
                v[i] += tdelta * a[i];

            auto dt_dtau { 1.0 / sqrt(-k.magnitude(v, x)) };
            vec v4; //four-velocity (v4[0] = dt / dtau
            for(size_t i{}; i < 4; ++i) v4[i] = v[i] * dt_dtau;
            auto vm { k.magnitude(v4, x) };

            if(abs(vm + 1) > 1e-3)
            {
                printf("Numeric instability too high!\n");
                break;
            }

            if (cnt % 100 == 0)
            {
                printf("X = [%4.4lf, %4.4lf, %4.4lf, %4.4lf]\t", x[0], x[1], x[2], x[3]);
                printf("dX/dτ = [%4.4lf, %4.4lf, %4.4lf, %4.4lf] mag = %4.4lf\n", v4[0], v4[1], v4[2], v4[3], vm);
            }

            ++cnt;
        }
        /*for (size_t i{}; i < 4; ++i)
            for (size_t j{}; j < 4; ++j)
                for (size_t k{}; k < 4; ++k)
                    std::cout << GV[i][j][k] << ", ";*/
        std::cout << '\n';
    }
#endif

#ifdef POLAR
    void polar()
    {
        using namespace symb;
        auto r{ v<'r'> };
        auto a{ v<'a'> };

        auto vars{ arr(r, a) };

        auto map{ arr(r * Cos(a), r * Sin(a)) };

        constexpr auto D{ vars.size() };

        //Jacobian
        auto J{ For<'v', D>(For<'u', D>(d(Ix<'v'>(map)) / d(Ix<'u'>(vars)))) };

        auto g{ For<'a', D>(For<'b', D>(Sum<'c', D>((Ix<'a'>(Ix<'c'>(J)) * Ix<'b'>(Ix<'c'>(J)))))) };
        auto G{ Christoffel(g, vars) };
        auto R{ Riemann(G, vars) };

        std::cout << "Metric: \n" << to_str(g.to_str()) << '\n';
        std::cout << "Γ:      \n" << to_str(G.to_str()) << '\n';
        std::cout << "R:      \n" << to_str(R.to_str()) << '\n';
    }
#endif

#ifdef SPHERE_SURFACE
    void sphere_surface()
    {
        using namespace symb;
        auto nul{ c<0.0> };

        auto a{ v<'a'> }; //longitude
        auto b{ v<'b'> }; //lattitude
        auto r{ v<'r'> };

        auto vars{ arr(a, b) };

        auto g{ arr(arr(r * r, nul), arr(nul, r * r * Sin(a) * Sin(a))) };
        auto G{ Christoffel(g, vars) };
        auto R{ Riemann(G, vars) };

        std::cout << "Metric: \n" << to_str(g.to_str()) << '\n';
        std::cout << "Γ:      \n" << to_str(G.to_str()) << '\n';
        std::cout << "R:      \n" << to_str(R.to_str()) << '\n';
    }
#endif

#ifdef MATRIX_INVERSE
    void matrix_inverse_3()
    {  
        using namespace symb;
        auto matrix{ arr(arr(v<'a'>, v<'b'>, v<'c'>), arr(v<'d'>, v<'e'>, v<'f'>), arr(v<'g'>, v<'h'>, v<'i'>)) };
        auto inverse{ matrix_inverse(matrix) };

        std::cout << to_str(matrix.to_str()) << '\n';
        std::cout << to_str(inverse.to_str()) << '\n';
    }                
#endif

#ifdef MATRIX_INVERSE_4
    void matrix_inverse_4()
    {  
        using namespace symb;
        auto matrix{ arr(arr(v<'a'>, v<'b'>, v<'c'>, v<'d'>), arr(v<'e'>, v<'f'>, v<'g'>, v<'h'>), arr(v<'i'>, v<'j'>, v<'k'>, v<'l'>), arr(v<'m'>, v<'n'>, v<'o'>, v<'p'>)) };
        auto inverse{ matrix_inverse(matrix) };

        std::cout << to_str(matrix.to_str()) << '\n';
        std::cout << to_str(inverse.to_str()) << '\n';
    }                
#endif
}


int main()
{
#ifdef WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

#ifdef MATRIX_INVERSE
    matrix_inverse_3();
#endif

#ifdef MATRIX_INVERSE_4
    matrix_inverse_4();
#endif

#ifdef SCHWARZSCHILD
    schwarzschild();
#endif

#ifdef POLAR
    polar();
#endif

#ifdef SPHERE_SURFACE
    sphere_surface();
#endif

#ifdef KERR
    kerr_geodesic();
#endif

#if 0
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
#endif
}
