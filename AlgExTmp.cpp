#include <iostream>
#include <locale>
#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif
#include <array>
#include "symb.h"
#include "matrix.h"
#include "riemann.h"
#include "kerr.h"
#include "str.h"
#include "vis.h"
#include "affine_space.h"
#include "riemann.h"

#include <fstream>
#include <random>

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

    void kerr_field()
    {
        using vec = affine_space<4>::vec;
        using point = affine_space<4>::point;

        struct mark
        {
            point p;
            vec v{ 1.0, 0.0, 0.0, 0.0 };
        };

        vis::init();

        double J{ 1 };
        double M{ 1 };
        kerr k { M, J };

        std::vector<mark> ms;
        std::vector<vis::point> glps;

        for (size_t xi{}; xi < 50; ++xi)
            for (size_t yi{}; yi < 50; ++yi)
            {
                auto xf{ (xi - 25.0) / 1.0 };
                auto yf{ (yi - 25.0) / 1.0 };
                auto r{ sqrt(xf * xf + yf * yf) };
                auto phi { atan2(xf, yf) };

                ms.emplace_back(point { 0, r, 3.14159265359 / 2, phi });
            }

        double delta{ 5e-1 };
        for (size_t i{}; i < 100; ++i)
        {
            for (auto& m : ms)
            {
                m.v = m.v + riemann<4>::coordinate_acceleration(k.christoffel(m.p), m.v) * delta;
                m.p = m.p + m.v * delta;
            }
        }

        std::transform(ms.begin(), ms.end(), std::back_inserter(glps),
            [&k](const mark& m) {
                return vis::point{ k.sphere_to_cart(m.p).coords, { 1., 0., 0., 1. } };
            }
        );
        vis::draw_points(glps.data(), glps.size());
        vis::swap();
        fgetchar();
    }

#if 0
    void kerr_geodesic()
    {
        vis::init();

        double J{ 0 };
        double M{ 1 };
        kerr k { M, J };

        std::minstd_rand rnd;
        std::normal_distribution dist_r { 40., 1e-1 };
        std::normal_distribution dist_r2 { 9., 1e-1 };
        std::uniform_real_distribution dist_p { 0., 3.14159265359 * 2 };
        std::uniform_real_distribution dist_t { 0., 3.14159265359 };
        std::normal_distribution dist_vr { 0., 3e-5 };
        std::normal_distribution dist_vt { 0., 0. };
        std::normal_distribution dist_vp { -9e-3, 9e-4 };

        std::vector<ppoint> ps;
        std::vector<vis::point> glps;
        for (size_t i{}; i < 1000; ++i)
        {
#if 1
            ps.emplace_back(
                vec { 0, dist_r2(rnd), 3.14159265359 / 2, dist_p(rnd) },
                vec { 1, 0, 0, 0.01 }
            );
#else
            ps.emplace_back(
                vec { 0, dist_r(rnd), dist_t(rnd), dist_p(rnd) },
                vec { 1.0, dist_vr(rnd), dist_vt(rnd), dist_vp(rnd) }
            );
#endif
            glps.emplace_back();
        }
//        std::array<double, 4> x { {0., 20, 3.1415926536 / 2 + 1e-2, 0 } };
//        std::array<double, 4> v { {1.0, 0., 0., 3e-3} }; // coordinate velocity (v[0] = dt/dt = 1)

        size_t cnt{};
        for (;;)
        {
            double tdelta{ 1e-1 };
            for (size_t i{}; i < ps.size(); ++i)
            {
                auto& p{ ps[i] };
                auto& glp{ glps[i] };

                auto GV{ k.christoffel(p.x) };
                auto a{ riemann<4>::coordinate_acceleration(GV, p.v) };

                for (size_t i{}; i < 4; ++i)
                    p.x[i] += tdelta * p.v[i];

                for (size_t i{}; i < 4; ++i)
                    p.v[i] += tdelta * a[i];

                auto dtau_dt{ sqrt(-k.magnitude(p.v, p.x)) };
                auto dt_dtau{ 1.0 / dtau_dt };
                vec v4; //four-velocity (v4[0] = dt / dtau
                for (size_t i{}; i < 4; ++i) v4[i] = p.v[i] * dt_dtau;

                auto vm{ k.magnitude(v4, p.x) };

                glp.pos = k.cyl_to_cart(p.x);
                if (isnan(dtau_dt))
                    dtau_dt = 0;
                auto clamped_dtau_dt{ std::clamp(dtau_dt, 0., 1.) };
                glp.color = { clamped_dtau_dt, 0, 1 - clamped_dtau_dt, 1 };

             /*   auto& x{ p.x };
                auto& xc{ glp.pos };
                auto& v{ p.v };
                printf("\33[2K\rT = %4.4lf X = Cyl[%4.4lf, %4.4lf, %4.4lf] Cart[%4.4lf, %4.4lf, %4.4lf] ", x[0], x[1], x[2], x[3], xc[1], xc[2], xc[3]);
                printf("dτ/dt = %4.4lf dX/dτ = [%4.4lf, %4.4lf, %4.4lf, %4.4lf] mag = %4.4lf", dtau_dt, v4[0], v4[1], v4[2], v4[3], vm);
                fflush(stdout);*/
                
            }
            vis::draw_points(glps.data(), glps.size());
            vis::swap();

            ++cnt;
        }
    }
#endif
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
        auto V{ v<'V'> }; //longitude
        auto r{ v<'r'> };

        auto vars{ arr(a, b) };

        auto g{ arr(arr(r * r, nul), arr(nul, r * r * Sin(a) * Sin(a))) };
        auto G{ Christoffel(g, vars) };
        auto R{ Riemann(G, vars) };
        auto ds{ Sqrt(Sum<'u', 2>(Sum<'v', 2>(Ix<'v'>(Ix<'u'>(g)) * Ix<'u'>(V) * Ix<'v'>(V)))) };

        std::cout << "Metric: \n" << to_str(g.to_str()) << '\n';
        std::cout << "Γ:      \n" << to_str(G.to_str()) << '\n';
        std::cout << "R:      \n" << to_str(R.to_str()) << '\n';

        using vec2 = std::array<double, 2>;
        double radius{ 1.0 };
        vec2 p{ 1.1, 0 }, p_{ 0, 1 }, p__{};
        using tens2_3 = std::array<std::array<std::array<double, 2>, 2>, 2>;
        tens2_3 GV;
        double s{};

        auto ctx{
            [&](char16_t v, size_t ix) {
                switch (v) {
                case 'r': return radius;
                case 'a': return p[0];
                case 'b': return p[1];
                case 'V': return p_[ix];
                }
            }
        };

        for (size_t cnt{};;++cnt)
        {
            auto tdelta{ 1e-8 };
            G.eval_array(ctx, &GV[0][0][0]);

            std::fill(p__.begin(), p__.end(), 0);
            for (int i{}; i < 2; ++i)
                for (int j{}; j < 2; ++j)
                    for (int k{}; k < 2; ++k)
                        p__[i] += -GV[i][j][k] * p_[j] * p_[k];

            for (std::size_t i{}; i < 2; ++i)
                p_[i] += tdelta * p__[i];

            for (std::size_t i{}; i < 2; ++i)
                p[i] += tdelta * p_[i];

            s += tdelta * ds.eval(ctx);

/*
            std::cout << "{";
            for (int i{}; i < 2; ++i)
            {
                if (i) std::cout << ", ";
                std::cout << "{";
                for (int j{}; j < 2; ++j)
                {
                    if (j) std::cout << ", ";
                    std::cout << "{";
                    for (int k{}; k < 2; ++k)
                    {
                        if (k) std::cout << ", ";
                        std::cout << GV[i][j][k];
                    }
                    std::cout << "}";
                }
                std::cout << "}";
            }
            std::cout << "}\n";
*/
            if (cnt % 1000000 == 0)
            {
                printf("\33[2K\rX = [%4.4lf, %4.4lf] V = [%4.4lf, %4.4lf] S = %4.4lf", p[0], p[1], p_[0], p_[1], s);
                fflush(stdout);
               // fgetchar();
            }
        }
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
    kerr_field();
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
