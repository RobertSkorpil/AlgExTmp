#include "kerr.h"
#ifndef __INTELLISENSE__
#include "symb.h"
#endif
#include "riemann.h"
#include <iostream>

using namespace symb;
namespace
{
    struct math_t
    {
      constexpr static auto t{ v<'t'> };
      constexpr static auto r{ v<'r'> };
      constexpr static auto theta{ v<u'Θ'> };
      constexpr static auto phi{ v<u'φ'> };

      constexpr static auto J{ v<'J'> };
      constexpr static auto M{ v<'M'> };
      constexpr static auto _0{ c<0.0> };
      constexpr static auto _1{ c<1.0> };
      constexpr static auto _2{ c<2.0> };
      constexpr static auto a{ J / (_2 * M) };
      constexpr static auto S{ r * r + a * a * Cos(theta) * Cos(theta) };
      constexpr static auto D{ r * r - _2 * M * r + a * a };
      constexpr static auto g{ arr(
          arr(-(_1 - _2 * M * r / S), _0, _0, -_2 * M * r * a * Sin(theta) * Sin(theta) / S),
          arr(_0, S / D, _0, _0),
          arr(_0, _0, S, _0),
          arr(-_2 * M * r * a * Sin(theta) * Sin(theta) / S, _0, _0, (r * r + a * a + _2 * M * r * a * a * Sin(theta) * Sin(theta) / S) * Sin(theta) * Sin(theta))
      ) };
      
      constexpr static auto A { v<'A'> };
      constexpr static auto magn { Sum<'i', 4>(Sum<'j', 4>(Ix<'j'>(Ix<'i'>(g)) * Ix<'i'>(A) * Ix<'j'>(A))) };

      constexpr static auto vars{ arr(t, r, theta, phi) };
      constexpr static auto G{ Christoffel(g, vars) };
    };
    
    struct context
    {
      double J, M, t, r, theta, phi;
      vec A;

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
          case 'A':
              return A[*index];
          default:
              return std::numeric_limits<double>::quiet_NaN();
          }
      }
    };

    std::function<double(vec v, point p)>
    make_magnitude(double M_, double J_)
    {
      return [=](vec v, point p)
      {
        context ctx{ J_, M_, p[0], p[1], p[2], p[3], v };
        return math_t::magn.eval(ctx);
      };
    }


    std::function<tensor3(point p)>
    make_christoffel(double M_, double J_)
    {

      std::cout << "Metric: \n" << to_str(math_t::g.to_str()) << '\n';
      std::cout << "Γ:      \n" << to_str(math_t::G.to_str()) << '\n';

      return [=](point p)
      {
          context ctx{ J_, M_, p[0], p[1], p[2], p[3] };
          tensor3 val;
          math_t::G.eval_array(ctx, &val[0][0][0]);
          return val;
      };
    }
}

kerr::kerr(double M, double J)
  : christoffel { make_christoffel(M, J) }, magnitude { make_magnitude(M, J) } {}
