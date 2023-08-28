#pragma once
#include <type_traits>
#include <string>
#include <optional>
#include <algorithm>
#include <math.h>
#include "str.h"


#ifdef _MSC_VER
#define INLINE __forceinline
#else
#define INLINE
#endif

using scalar = double;

template<size_t N>
struct StringLiteral {
    constexpr StringLiteral(const char(&str)[N]) {
        std::copy_n(str, N - 1, value);
    }

    char value[N - 1];

    constexpr std::u16string to_str() const
    {
        return std::u16string{ std::begin(value), std::end(value) };
    }

    constexpr auto hash() const
    {
        size_t v{};
        for (size_t i{}; i < N - 1; ++i)
            v += value[i];
        return v;
    }
};

namespace symb
{
    struct index
    {
        char16_t n{};

        constexpr index() = default;
        constexpr index(const index&) = default;
        constexpr index(char16_t n) : n{ n } {}
        constexpr bool operator ==(const index&) const = default;
        constexpr operator bool() const { return n; }
    };

    struct var
    {
        char16_t n;
        size_t ix;

        constexpr var(const var&) = default;
        constexpr var(char16_t n) : n{ n }, ix { std::numeric_limits<size_t>::max() } {}
        constexpr var(char16_t n, size_t ix) : n{ n }, ix{ ix } {}

        constexpr bool operator ==(const var&) const = default;
    };

    struct dummy_eval_context
    {
        auto var_value(const var& v, std::optional<size_t> index) const
        {
            return 0.0;
        }

        void set_var_value(const var& v, std::optional<size_t> index, scalar value) const
        {
        }
    };

    struct no_index {
        using terminator = int;
    };

    template<typename T>
    concept ActualIndexAssignment = requires { T::ix; T::val; };

    template<typename T>
    concept IndexAssignment = ActualIndexAssignment<T> || std::is_same_v<T, no_index>;

    template<typename T>
    concept VarExp = requires { std::remove_reference_t<T>::is_var; };

    template<typename T>
    concept Expr = requires { std::remove_reference_t<T>::eval(dummy_eval_context{}); };

    template<typename T>
    concept Valuation = requires { T::var_value('x'); };

    template<typename T>
    concept ConstantLike = requires { T::is_const; }&& Expr<T>;

    template<typename T>
    concept NonConstant = requires { T::non_const; }&& Expr<T>;

    template<typename T>
    concept NonPower = requires { T::non_power; }&& NonConstant<T>;

    template<typename T>
    concept ActualAction = requires { T::is_action; };

    template<typename T>
    concept Action = ActualAction<T> || Expr<T>;

    template<typename T>
    concept DifferentialLike = requires { T::is_differential; };

    template<Expr ExprT>
    constexpr auto expand(ExprT expr);

    template<Expr ExprT>
    auto simplify(ExprT expr);

    template<Expr ExprT>
    constexpr size_t hash = 0;

    template<index Ix, size_t Val, IndexAssignment Next>
    struct index_assignment
    {
        static constexpr index ix { Ix };
        static constexpr size_t val { Val };
        using next = Next;
    };
    
    template<index Ix, size_t Val, IndexAssignment Old>
    constexpr auto add_index_assignment()
    {
        if constexpr (std::is_same_v<Old, no_index>)
            return index_assignment<Ix, Val, Old>{};
        else if constexpr (Ix.n < Old::ix.n)
            return index_assignment<Ix, Val, Old>{};
        else if constexpr (std::true_type{})
            return index_assignment<Old::ix, Old::val, decltype(add_index_assignment<Ix, Val, typename Old::next>())>{};
    }

    template<Expr ExprT>
    constexpr bool has_indexed_expr();

    template<Expr ExprT>
    constexpr bool all_index_assigned();

    template<IndexAssignment Ix, index ix>
    constexpr bool is_index_assigned();

    template<IndexAssignment Ix, index ix>
    constexpr bool is_index_assigned()
    {
        if constexpr (std::is_same_v<Ix, no_index>)
            return false;
        else if constexpr (Ix::ix == ix)
            return true;
        else
            return is_index_assigned<typename Ix::next, ix>();
    }

    template<IndexAssignment Ix, index ix>
    constexpr size_t ix_val_static();

    template<IndexAssignment Ix, index ix>
    constexpr size_t ix_val_static()
    {
        if constexpr (std::is_same_v<Ix, no_index>)
            return ~0ull;
        else if (Ix::ix == ix)
            return Ix::val;
        else
            return ix_val_static<typename Ix::next, ix>();
    }


    template<IndexAssignment Ix>
    constexpr std::optional<size_t> ix_val(index ix);

    template<IndexAssignment Ix>
    constexpr std::optional<size_t> ix_val(index ix)
    {
        if constexpr (std::is_same_v<Ix, no_index>)
            return {};
        else if (Ix::ix == ix)
            return Ix::val;
        else
            return ix_val<typename Ix::next>(ix);
    }

    template<Expr ExprT, index ix, IndexAssignment Ix = no_index>
    struct IndexedExpr
    {
        static constexpr auto non_power{ true };
        using expr_t = ExprT;        
        using ix_assign = Ix;
        static constexpr index expr_ix = ix;

        template<typename Valuation>
        static auto INLINE eval(const Valuation& ctx)
        {
            return std::numeric_limits<scalar>::quiet_NaN();
        }

        template<VarExp D>
        static constexpr auto diff()
        {
            return ExprT::diff();
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        template<Expr RHS>
        constexpr auto operator =(const RHS& rhs);

        static std::u16string to_str()
        {
            return ExprT::to_str() + u"[" + ix.n + u"]";
        }
    };

    template<typename T>
    struct is_indexed_expr : std::false_type {};

    template<Expr ExprT, index ix, IndexAssignment Ix>
    struct is_indexed_expr<IndexedExpr<ExprT, ix, Ix>> : std::true_type {};

    using func_single_t = scalar(scalar);

    template<scalar c, IndexAssignment Ix = no_index>
    struct Constant
    {
        static constexpr auto non_power{ true };
        static constexpr auto value = c;
        using ix_assign = Ix;
        template<typename Valuation>
        static auto INLINE eval(const Valuation& ctx)
        {
            return c;
        }

        template<VarExp D>
        static auto diff()
        {
            return Constant<0.0>{};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            std::string s { abs(c - floor(c)) < 1e-5 ? std::to_string((int)c) : std::to_string(c) };
            return to_u16str(s);
        }
    };

    using ZeroExpr = Constant<0.0>;
    using OneExpr = Constant<1.0>;

    template<typename T>
    struct is_constant : std::false_type {};

    template<scalar c, IndexAssignment Ix>
    struct is_constant<Constant<c, Ix>> : std::true_type {};

    template<scalar c, IndexAssignment Ix>
    struct is_constant<const Constant<c, Ix>> : std::true_type {};

    template<var vr, IndexAssignment Ix = no_index>
    struct VarExpr
    {
        static constexpr auto non_power{ true };
        static constexpr bool is_var{ true };
        static constexpr bool non_const{ true };
        static constexpr var V { vr };
        using ix_assign = Ix;

        template<typename Valuation>
        static auto INLINE eval(const Valuation& ctx)
        {
            return ctx.var_value(vr, vr.ix);
        }

        template<VarExp D>
        constexpr static auto diff()
        {
            constexpr auto d{ D::V };
            if constexpr (d == vr)
                return OneExpr{};
            else if constexpr (d.n == vr.n)
            {
                constexpr auto ix1{ ix_val<Ix>(vr.ix) };
                constexpr auto ix2{ ix_val<typename D::ix_assign>(D::V.ix) };
                if constexpr (ix1 && ix2 && *ix1 == *ix2)
                    return OneExpr{};
                else
                    return ZeroExpr{};
            }
            else if constexpr (d != vr)
                return ZeroExpr{};
            else
                return ZeroExpr{};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            return std::u16string{ vr.n } + ((vr.ix != std::numeric_limits<size_t>::max()) ? std::u16string{ wchar_t(u'₀' + vr.ix) } : u"");
        }

        template<Expr RHS>
        constexpr auto operator =(const RHS& rhs);

        static constexpr size_t hash()
        {
            return V.n;
        }
    };

    template<var vr, IndexAssignment Ix>
    constexpr size_t hash<VarExpr<vr, Ix>> = VarExpr<vr, Ix>::hash();

    template<typename T>
    struct is_var_expr : std::false_type {};

    template<var vr, IndexAssignment Ix>
    struct is_var_expr<VarExpr<vr, Ix>> : std::true_type{};

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    struct Assignment
    {
        static constexpr auto non_power{ true };
        static constexpr bool is_action{ true };
        using ix_assign = Ix;

        template<typename Valuation>
        static void INLINE eval(const Valuation& ctx)
        {
            ctx.set_var_value(LHS::V, ix_val<LHS::Ix>(LHS::V.ix), RHS::eval(ctx));
        }

        static std::u16string to_str()
        {
            return std::remove_reference_t<LHS>::to_str() + u" = " + RHS::to_str();
        }
    };

    template<Expr ExprT, index ix, IndexAssignment Ix>
    template<Expr RHS> 
    constexpr auto IndexedExpr<ExprT, ix, Ix>::operator =(const RHS &rhs)
    {
        return Assignment<decltype(*this), RHS, Ix>{};
    }

    struct nothing
    {
        using act_t = nothing;
        using next_t = nothing;

        static constexpr bool is_action{ true };

        template<typename Valuation>
        static void INLINE eval(const Valuation& ctx)
        {}

        template<typename Valuation>
        static double * eval_array(const Valuation& ctx, double * out)
        {
            return out;
        }

        static std::u16string to_str()
        {
            return {};
        }

        static std::u16string to_str_internal()
        {
            return {};
        }

        static constexpr size_t size()
        {
            return 0;
        }
    };


    template<typename T>
    struct is_array : std::false_type {};

    template<Action Act, Action Next = nothing, IndexAssignment Ix = no_index>
    struct Array;

    template<Action Act, Action Next, IndexAssignment Ix>
    struct is_array<Array<Act, Next, Ix>> : std::true_type{};

    template<Action Act, Action Next, IndexAssignment Ix>
    struct Array
    {
        static constexpr auto non_power{ true };
        static constexpr int non_const{};
        static constexpr bool is_action{ true };
        using act_t = Act;
        using next_t = Next;
        using ix_assign = Ix;

        template<typename Valuation>
        static void INLINE eval(const Valuation& ctx)
        {
            Act::eval(ctx);
            Next::eval(ctx);
        }

        template<typename Valuation>
        static double * eval_array(const Valuation& ctx, double *out)
        {
            if constexpr (is_array<Act>{})
                out = Act::eval_array(ctx, out);
            else
                *out++ = Act::eval(ctx);

            out = Next::eval_array(ctx, out);
            return out;
        }

        template<VarExp D>
        static constexpr auto diff()
        {
            return Array{};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str_internal()
        {
            return Act::to_str() + (std::is_same_v<Next, nothing> ? u"" :  u"," + Next::to_str_internal());
        }

        static std::u16string to_str()
        {
            return u"{" + to_str_internal() + u"}";
        }

        static constexpr size_t size()
        {
            return 1 + Next::size();
        }
    };

    template<var v, IndexAssignment Ix>
    template<Expr RHS>
    constexpr auto VarExpr<v, Ix>::operator =(const RHS& rhs)
    {
        return Assignment<decltype(*this), RHS, Ix>{};
    }

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    struct SumExpr;

    template<typename T>
    struct is_sum_expr : std::false_type {};

    template<Expr LHS, Expr RHS, IndexAssignment Ix>
    struct is_sum_expr<SumExpr<LHS, RHS, Ix>> : std::true_type {};

    template<Expr LHS, Expr RHS, IndexAssignment Ix>
    struct SumExpr
    {
        static constexpr auto non_power{ true };
        static constexpr int non_const{};
        using lhs_t = LHS;
        using rhs_t = RHS;
        using ix_assign = Ix;

        template<typename Valuation>
        static auto INLINE eval(const Valuation& ctx)
        {
            return LHS::eval(ctx) + RHS::eval(ctx);
        }

        template<VarExp D>
        static constexpr auto diff()
        {
            return SumExpr<typename LHS::template diff_t<D>, typename RHS::template diff_t<D>, Ix>{};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            return u"(" + LHS::to_str() + u" + " + RHS::to_str() + u")";
        }
    };

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    struct ProdExpr
    {
        static constexpr int non_const{};
        using lhs_t = LHS;
        using rhs_t = RHS;
        using ix_assign = Ix;

        template<typename Valuation>
        static auto INLINE eval(const Valuation& ctx)
        {
            return LHS::eval(ctx) * RHS::eval(ctx);
        }

        template<VarExp D>
        static constexpr auto diff()
        {
            return SumExpr<ProdExpr<LHS, typename RHS::template diff_t<D>>, ProdExpr<typename LHS::template diff_t<D>, RHS>>{};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            return u"(" + LHS::to_str() + u" * " + RHS::to_str() + u")";
        }
    };

    template<Expr LHS, Expr RHS, IndexAssignment Ix>
    constexpr size_t hash<ProdExpr<LHS, RHS, Ix>> = std::numeric_limits<size_t>::max();

    template<typename T>
    struct is_prod_expr : std::false_type {};

    template<Expr LHS, Expr RHS, IndexAssignment Ix>
    struct is_prod_expr<ProdExpr<LHS, RHS, Ix>> : std::true_type {};

    template<typename f, Expr ArgExpr, IndexAssignment Ix = no_index>
    struct FuncExpr;

    template<func_single_t f, StringLiteral name>
    struct Func
    {
        static double eval(double arg)
        {
            return f(arg);
        }

        template<Expr ArgExpr>
        constexpr auto operator()(const ArgExpr& arg)
        {
            return FuncExpr<Func, ArgExpr>{};
        }

        static std::u16string to_str()
        {
            return name.to_str();
        }

        static constexpr size_t hash()
        {
            return name.hash();
        }
    };

    template<Expr ExprT>
    struct Differential
    {
        static constexpr bool is_differential{ true };
        using expr_t = ExprT;
    };

    template<Expr ExprT>
    constexpr auto d(ExprT expr)
    {
        return Differential<ExprT>{};
    }

    template<typename T> struct test;

    template<Expr F, Expr V, IndexAssignment Ix = no_index>
    struct DerivativeExpr
    {
        static constexpr auto non_power{ true };
        static constexpr int non_const{};
        using f_t = F;
        using v_t = V;
        using ix_assign = Ix;

        template<typename Valuation>
        static auto INLINE eval(const Valuation& ctx)
        {
            return std::numeric_limits<scalar>::quiet_NaN();
        }

        constexpr static auto derive()
        {
            if constexpr (!is_indexed_expr<F>{} && is_var_expr<V>{})
                return typename F::template diff_t<V>{};
            else
                return DerivativeExpr{};
        }

        static std::u16string to_str()
        {
            return u" d(" + F::to_str() + u")/d(" + V::to_str() + u") ";
        }
    };

    template<typename T>
    struct is_derivative_expr : std::false_type {};

    template<Expr F, Expr V, IndexAssignment Ix>
    struct is_derivative_expr<DerivativeExpr<F, V, Ix>> : std::true_type {};

    template<DifferentialLike F, DifferentialLike V>
    constexpr auto operator /(const F&, const V&)
    {
        using DE = DerivativeExpr<typename F::expr_t, typename V::expr_t>;
        return DE::derive();
    }

    template<typename f, Expr ArgExpr, IndexAssignment Ix>
    struct func_derivative {};

    using sin_t = Func<sin, "Sin">; 
    using cos_t = Func<cos, "Cos">;

    scalar do_inverse(scalar x);
    using inv_t = Func<do_inverse, "inv">;

    extern sin_t Sin;
    extern cos_t Cos;
    extern inv_t Inv;

    template<Expr ArgExpr, IndexAssignment Ix>
    struct func_derivative<sin_t, ArgExpr, Ix>
    {
        using expr_t = FuncExpr<cos_t, ArgExpr, Ix>;
    };

    template<Expr ArgExpr, IndexAssignment Ix>
    struct func_derivative<cos_t, ArgExpr, Ix>
    {
        using expr_t = ProdExpr<Constant<-1.0>, FuncExpr<sin_t, ArgExpr, Ix>, Ix>;
    };

    template<Expr ArgExpr, IndexAssignment Ix>
    struct func_derivative<inv_t, ArgExpr, Ix>
    {
        using expr_t = FuncExpr<inv_t, ProdExpr<ProdExpr<Constant<-1.0>, ArgExpr, Ix>, ArgExpr, Ix>, Ix>;
    };

    template<typename F, Expr ArgExpr, IndexAssignment Ix>
    struct FuncExpr
    {
        static constexpr auto non_power{ true };
        static constexpr int non_const{};
        using f_t = F;
        using arg_t = ArgExpr;
        using ix_assign = Ix;

        template<typename Valuation>
        static auto INLINE eval(const Valuation& ctx)
        {
            return F::eval(ArgExpr::eval(ctx));
        }

        template<VarExp D>
        static constexpr auto diff()
        {
            return ProdExpr< typename func_derivative<F, ArgExpr, Ix>::expr_t, typename ArgExpr::template diff_t<D>, Ix>{};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            return F::to_str() + u"[" + ArgExpr::to_str() + u"]";
        }

        static constexpr size_t hash()
        {
            return F::hash();
        }
    };

    template<Expr ExprT, int Power, IndexAssignment Ix = no_index>
    struct PowerExpr
    {
        static constexpr int non_const{};
        static constexpr auto power = Power;
        using expr_t = ExprT;
        using ix_assign = Ix;

        template<typename Valuation>
        static auto INLINE eval(const Valuation& ctx)
        {
            if constexpr (Power == 0)
                return 1.0;
            else if constexpr (Power > 0)
                return ExprT::eval(ctx) * PowerExpr<ExprT, Power - 1, Ix>::eval(ctx);
            else if constexpr (Power < 0)
                return 1.0 / ExprT::eval(ctx) * PowerExpr<ExprT, Power + 1, Ix>::eval(ctx);
        }

        template<VarExp D>
        static constexpr auto diff()
        {
            return to_product().template diff<D>();
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            return ExprT::to_str() + u"^" + to_u16str(std::to_string(Power));
        }

        constexpr static auto to_product_abs()
        {
            if constexpr (Power == 0)
                return OneExpr{};
            else if constexpr (std::true_type{})
                return ProdExpr<ExprT, decltype(PowerExpr<ExprT, Power - 1, Ix>::to_product_abs()), Ix>{};
        }

        constexpr static auto to_product()
        {
            if constexpr (Power < 0)
                return FuncExpr<inv_t, decltype(PowerExpr<ExprT, -Power, Ix>::to_product_abs()), Ix>{};
            else if constexpr (std::true_type{})
                return to_product_abs();
        }
    };

    template<typename T>
    struct is_power_expr : std::false_type {};

    template<Expr ExprT, int Power, IndexAssignment Ix>
    struct is_power_expr<PowerExpr<ExprT, Power, Ix>> : std::true_type {};

    template<Expr ExprT, int power>
    constexpr size_t hash<PowerExpr<ExprT, power>> = hash<ExprT>;

    template<typename F, Expr ArgExpr, IndexAssignment Ix>
    constexpr size_t hash<FuncExpr<F, ArgExpr, Ix>> = FuncExpr<F, ArgExpr, Ix>::hash();

    template<size_t i, size_t v, Expr Main, Expr ExprT>
    constexpr auto do_index_seq(ExprT expr);

    template<index ix, Expr ExprT>
    constexpr auto Ix(ExprT expr)
    {
        return IndexedExpr<ExprT, ix>{};
    }

    template<size_t ix, Expr ExprT>
    constexpr auto I(ExprT expr)
    {
        return do_index_seq<0, ix, ExprT>(expr);
    }

    template<typename T>
    struct is_func_expr : std::false_type {};

    template<typename F, Expr ArgExpr, IndexAssignment Ix>
    struct is_func_expr<FuncExpr<F, ArgExpr, Ix>> : std::true_type {};

    template<Expr ExprT, IndexAssignment Ix = no_index>
    constexpr auto operator -(ExprT expr)
    {
        return ProdExpr<Constant<-1.0>, ExprT, Ix>{};
    }

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    constexpr auto operator -(LHS lhs, RHS rhs)
    {
        return SumExpr<LHS, ProdExpr<Constant<-1.0>, RHS, Ix>, Ix>{};
    }

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    constexpr auto operator +(LHS lhs, RHS rhs)
    {
        return SumExpr<LHS, RHS, Ix>{};
    }

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    constexpr auto operator *(LHS lhs, RHS rhs)
    {
        return ProdExpr<LHS, RHS, Ix>{};
    }

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    constexpr auto operator /(LHS lhs, RHS rhs)
    {
        return ProdExpr<LHS, PowerExpr<RHS, -1, Ix>, Ix>{};
    }

    template<scalar val>
    constexpr Constant<val> c{};

    template<char16_t vr, size_t ix>
    constexpr auto make_v()
    {
        constexpr var v{ vr, ix };
        return VarExpr<v>{};
    }

    template<char16_t vr, size_t ix = std::numeric_limits<size_t>::max()>
    constexpr auto v { make_v<vr, ix>() };

    constexpr nothing make_array() {
        return {};
    };

    template<Expr ExprT, Expr... Rest>
    constexpr auto make_array(ExprT e, Rest... args)
    {
        return Array<ExprT, decltype(make_array(std::forward<Rest>(args)...))>{};
    }

    template<Expr... Args>
    constexpr auto arr(Args... rest) {
        return make_array(std::forward<Args>(rest)...);
    };

    template<index ix, size_t val, Expr ExprT>
    struct assign_index_t
    {
    };

    template<index ix, size_t val, IndexAssignment Next, Expr ExprT>
    struct assign_index_t<ix, val, IndexedExpr<ExprT, ix, Next>>
    {
        using expr_t = IndexedExpr<ExprT, ix, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr ExprT, index expr_ix>
    struct assign_index_t<ix, val, IndexedExpr<ExprT, expr_ix, Next>>
    {
        using expr_t = IndexedExpr<typename assign_index_t<ix, val, ExprT>::expr_t, expr_ix, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, scalar c>
    struct assign_index_t<ix, val, Constant<c, Next>>
    {
        using expr_t = Constant<c, Next>;
    };

    template<index ix, size_t val>
    struct assign_index_t<ix, val, nothing>
    {
        using expr_t = nothing;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, SumExpr<LHS, RHS, Next>>
    {
        using expr_t = SumExpr<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, ProdExpr<LHS, RHS, Next>>
    {
        using expr_t = ProdExpr<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, DerivativeExpr<LHS, RHS, Next>>
    {
        using expr_t = DerivativeExpr<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, Assignment<LHS&, RHS, Next>>
    {
        using expr_t = Assignment<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, Assignment<LHS, RHS, Next>>
    {
        using expr_t = Assignment<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr E, int Power>
    struct assign_index_t<ix, val, PowerExpr<E, Power, Next>>
    {
        using expr_t = PowerExpr<typename assign_index_t<ix, val, E>::expr_t, Power, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, var v>
    struct assign_index_t<ix, val, VarExpr<v, Next>>
    {
        using expr_t = VarExpr<v, no_index>;
    };

    template<index ix, size_t val, IndexAssignment Next, Action Act, Action NextAct>
    struct assign_index_t<ix, val, Array<Act, NextAct, Next>>
    {
        using expr_t = Array<typename assign_index_t<ix, val, Act>::expr_t, typename assign_index_t<ix, val, NextAct>::expr_t, decltype(add_index_assignment<ix, val, Next>())>;
    };

    template<index ix, size_t val, IndexAssignment Next, typename F, Expr Arg>
    struct assign_index_t<ix, val, FuncExpr<F, Arg, Next>>
    {
        using expr_t = FuncExpr<F, typename assign_index_t<ix, val, Arg>::expr_t, decltype(add_index_assignment<ix, val, Next>())>;
    };


    template<Expr ExprT>
    struct unassign_index_t
    {
    };

    template<IndexAssignment Next, scalar c>
    struct unassign_index_t<Constant<c, Next>>
    {
        using expr_t = Constant<c, Next>;
    };

    template<IndexAssignment Next, scalar c>
    struct unassign_index_t<const Constant<c, Next>>
    {
        using expr_t = Constant<c, Next>;
    };

    template<>
    struct unassign_index_t<nothing>
    {
        using expr_t = nothing;
    };

    template<IndexAssignment Next, Expr LHS, Expr RHS>
    struct unassign_index_t<SumExpr<LHS, RHS, Next>>
    {
        using expr_t = SumExpr<typename unassign_index_t<LHS>::expr_t, typename unassign_index_t<RHS>::expr_t, no_index>;
    };

    template<IndexAssignment Next, Expr LHS, Expr RHS>
    struct unassign_index_t<ProdExpr<LHS, RHS, Next>>
    {
        using expr_t = ProdExpr<typename unassign_index_t<LHS>::expr_t, typename unassign_index_t<RHS>::expr_t, no_index>;
    };

    template<IndexAssignment Next, Expr LHS, Expr RHS>
    struct unassign_index_t<DerivativeExpr<LHS, RHS, Next>>
    {
        using expr_t = DerivativeExpr<typename unassign_index_t<LHS>::expr_t, typename unassign_index_t<RHS>::expr_t, no_index>;
    };

    template<IndexAssignment Next, Expr LHS, Expr RHS>
    struct unassign_index_t<Assignment<LHS&, RHS, Next>>
    {
        using expr_t = Assignment<typename unassign_index_t<LHS>::expr_t, typename unassign_index_t<RHS>::expr_t, no_index>;
    };

    template<IndexAssignment Next, Expr LHS, Expr RHS>
    struct unassign_index_t<Assignment<LHS, RHS, Next>>
    {
        using expr_t = Assignment<typename unassign_index_t<LHS>::expr_t, typename unassign_index_t<RHS>::expr_t, no_index>;
    };

    template<IndexAssignment Next, int Power, Expr E>
    struct unassign_index_t<PowerExpr<E, Power, Next>>
    {
        using expr_t = PowerExpr<E, Power>;
    };

    template<IndexAssignment Next, var v>
    struct unassign_index_t<VarExpr<v, Next>>
    {
        using expr_t = VarExpr<v, no_index>;
    };

    template<IndexAssignment Next, Action Act, Action NextAct>
    struct unassign_index_t<Array<Act, NextAct, Next>>
    {
        using expr_t = Array<typename unassign_index_t<Act>::expr_t, typename unassign_index_t<NextAct>::expr_t, no_index>;
    };

    template<IndexAssignment Next, typename F, Expr Arg>
    struct unassign_index_t<FuncExpr<F, Arg, Next>>
    {
        using expr_t = FuncExpr<F, typename unassign_index_t<Arg>::expr_t, no_index>;
    };

    template<Expr ExprT>
    constexpr auto expand(ExprT expr);

    template<index ix, size_t val, size_t range, Expr ExprT>
    constexpr auto sum_sumand(ExprT expr);

    template<index ix, size_t val, size_t range, Expr ExprT>
    constexpr auto sum_sumand(ExprT expr)
    {
        if constexpr (val == range)
            return ZeroExpr{};
        else if constexpr (val == range - 1)
            return typename assign_index_t<ix, val, decltype(expr)>::expr_t{};
        else if constexpr (std::true_type{})
            return SumExpr<typename assign_index_t<ix, val, decltype(expr)>::expr_t, decltype(sum_sumand<ix, val + 1, range>(expr)), decltype(add_index_assignment<ix, val, typename ExprT::ix_assign>())>{};
    }

    template<index ix, size_t range, Expr ExprT>
    constexpr auto Sum(ExprT expr)
    {
        return simplify(expand(sum_sumand<ix, 0, range>(simplify(expand(expr)))));
    }

    template<index ix, size_t val, size_t range, Expr ExprT>
    auto fac_factor(ExprT expr);

    template<index ix, size_t val, size_t range, Expr ExprT>
    auto fac_factor(ExprT expr)
    {
        if constexpr (val == range)
            return ZeroExpr{};
        else if constexpr (val == range - 1)
            return typename assign_index_t<ix, val, decltype(expr)>::expr_t{};
        else
            return ProdExpr<typename assign_index_t<ix, val, decltype(expr)>::expr_t, decltype(fac_factor<ix, val + 1, range, ExprT>(expr)), decltype(index_assignment<ix, val, typename ExprT::ix_assign>())>{};
    }

    template<index ix, size_t range, Expr ExprT>
    auto Prod(ExprT expr)
    {
        return simplify(expand(fac_factor<ix, 0, range>(simplify(expand(expr)))));
    }

    template<index ix, size_t val, size_t range, Expr ExprT>
    constexpr auto make_for(ExprT expr);

    template<index ix, size_t val, size_t range, Expr ExprT>
    constexpr auto make_for(ExprT expr)
    {
        if constexpr (val == range)
            return nothing{};
        else
            return Array<typename assign_index_t<ix, val, decltype(expr)>::expr_t, decltype(make_for<ix, val + 1, range, ExprT>(expr))> {};
    }

    template<index ix, size_t range, Expr ExprT>
    constexpr auto For(ExprT expr)
    {
        return simplify(expand(make_for<ix, 0, range>(simplify(expand(expr)))));
    }

    template<index ix, size_t val, size_t skip, size_t range, Expr ExprT>
    constexpr auto make_for_skip(ExprT expr);

    template<index ix, size_t val, size_t skip, size_t range, Expr ExprT>
    constexpr auto make_for_skip(ExprT expr)
    {
        if constexpr (val == range)
            return nothing{};
        else
        {
            constexpr auto next{ skip == val + 1 ? val + 2 : val + 1 };
            return Array<typename assign_index_t<ix, val, decltype(expr)>::expr_t, decltype(make_for_skip<ix, next, skip, range, ExprT>(expr))> {};
        }
    }

    template<index ix, size_t range, size_t skip, Expr ExprT>
    constexpr auto ForSkip(ExprT expr)
    {
        if constexpr (skip == 0)
            return simplify(expand(make_for<ix, 1, range>(simplify(expand(expr)))));
        else if constexpr (skip == range - 1)
            return simplify(expand(make_for<ix, 0, range - 1>(simplify(expand(expr)))));
        else
            return simplify(expand(make_for_skip<ix, 0, skip, range>(simplify(expand(expr)))));
    }

    template<typename T>
    struct debug;

    template<size_t v>
    struct debug_size_t;

    template<Expr ExprT>
    constexpr auto do_expand(ExprT expr);

    template<size_t i, size_t v, Expr Main, Expr ExprT>
    constexpr auto do_index_seq(ExprT expr)
    {
        if constexpr (std::true_type{})
        {
            if constexpr (i == v)
                return typename ExprT::act_t{};
            else
                return do_index_seq<i+1, v, Main>(typename ExprT::next_t{});
        }
    }

    template<IndexAssignment Ix, Expr ExprT>
    constexpr auto do_index(ExprT expr);

    template<IndexAssignment Ix, Expr ExprT>
    constexpr auto do_index(ExprT expr)
    {
        if constexpr (is_array<typename ExprT::expr_t>{})
        {
            constexpr auto v{ ix_val_static<typename ExprT::ix_assign, ExprT::expr_ix>() };
            return do_index_seq<size_t(0), v, ExprT>(do_expand(typename ExprT::expr_t{}));
        }
        else if constexpr (is_var_expr<typename ExprT::expr_t>{})
        {
            constexpr auto val{ ix_val_static<Ix, ExprT::expr_ix>() };
            return VarExpr<var{ ExprT::expr_t::V.n, val }, no_index> {};
        }
        else if constexpr (is_indexed_expr<typename ExprT::expr_t>{})
        {
            auto r{ do_expand(typename ExprT::expr_t{}) };
            if constexpr (!std::is_same_v<decltype(r), typename ExprT::expr_t>)
                return do_index<Ix>(IndexedExpr<decltype(r), ExprT::expr_ix, typename ExprT::ix_assign>{});
            else if constexpr (std::true_type{})
                return expr;
        }
        else if constexpr (std::true_type{})
            return expr;
    }

    template<Expr ExprT>
    constexpr auto do_expand(ExprT expr)
    {
        if constexpr (is_array<ExprT>{})
            return Array<decltype(do_expand(typename ExprT::act_t{})), decltype(do_expand(typename ExprT::next_t{})), typename ExprT::ix_assign > {};
        else if constexpr (is_indexed_expr<ExprT>{})
        {
            if constexpr (is_index_assigned<typename ExprT::ix_assign, ExprT::expr_ix>())
            {
                auto r{ do_index<typename ExprT::ix_assign>(expr) };
                if constexpr (std::is_same_v<decltype(r), ExprT>)
                    return expr;
                else if constexpr (std::true_type{})
                    return do_expand(r);
            }
            else if constexpr (std::true_type{})
                return expr;
        }
        else if constexpr (is_sum_expr<ExprT>{})
            return SumExpr<decltype(do_expand(typename ExprT::lhs_t{})), decltype(do_expand(typename ExprT::rhs_t{})), typename ExprT::ix_assign > {};
        else if constexpr (is_prod_expr<ExprT>{})
            return ProdExpr<decltype(do_expand(typename ExprT::lhs_t{})), decltype(do_expand(typename ExprT::rhs_t{})), typename ExprT::ix_assign > {};
        else if constexpr (is_derivative_expr<ExprT>{})
            return DerivativeExpr<decltype(do_expand(typename ExprT::f_t{})), decltype(do_expand(typename ExprT::v_t{})), typename ExprT::ix_assign > {}.derive();
        else if constexpr (is_func_expr<ExprT>{})
            return FuncExpr<decltype(typename ExprT::f_t{}), decltype(do_expand(typename ExprT::arg_t{})), typename ExprT::ix_assign > {};
        else if constexpr (std::true_type{})
            return expr;
    }

    template<Expr ExprT>
    constexpr auto expand(ExprT expr)
    {
        using R = decltype(do_expand(expr));
        if constexpr (std::is_same_v<ExprT, R>)
            return expr;
        else
            return expand(R{});
    }

    template<Expr ExprT>
    struct simplify_t { using expr_t = ExprT; };

    template<>
    struct simplify_t<Constant<-0.0>>
    {
        using expr_t = ZeroExpr;
    };

    template<Expr ExprT>
    struct simplify_t<PowerExpr<ExprT, 0>>
    {
        using expr_t = OneExpr;
    };

    template<Expr ExprT>
    struct simplify_t<PowerExpr<ExprT, 1>>
    {
        using expr_t = ExprT;
    };

    template<Expr ExprT, int Power>
    struct simplify_t<PowerExpr<ExprT, Power>>
    {
        using expr_t = PowerExpr<typename simplify_t<ExprT>::expr_t, Power>;
    };

    template<scalar s, NonConstant Sumand, IndexAssignment Ix>
    struct simplify_t<SumExpr<Sumand, ProdExpr<Constant<s>, Sumand, Ix>, Ix>>
    {
        using expr_t = ProdExpr<Constant<s + 1>, typename simplify_t<Sumand>::expr_t, Ix>;
    };

    template<scalar s, NonConstant Sumand, IndexAssignment Ix>
    struct simplify_t<SumExpr<ProdExpr<Constant<s>, Sumand, Ix>, Sumand, Ix>>
    {
        using expr_t = ProdExpr<Constant<s + 1>, typename simplify_t<Sumand>::expr_t, Ix>;
    };

    template<scalar s, scalar s2, NonConstant Sumand>
    struct simplify_t<SumExpr<ProdExpr<Constant<s>, Sumand>, ProdExpr<Constant<s2>, Sumand>>>
    {
        using expr_t = ProdExpr<Constant<s + s2>, typename simplify_t<Sumand>::expr_t>;
    };

    template<scalar a, scalar b, NonConstant RRHS>
    struct simplify_t<ProdExpr<Constant<a>, ProdExpr<Constant<b>, RRHS>/*, std::enable_if_t<a - b != 0, no_index*/>>
    {
        using expr_t = ProdExpr<Constant<a * b>, RRHS>;
    };

    template<NonConstant LHS, scalar a, NonConstant RRHS>
    struct simplify_t<ProdExpr<LHS, ProdExpr<Constant<a>, RRHS>>>
    {
        using expr_t = ProdExpr<Constant<a>, ProdExpr<LHS, RRHS>>;
    };

    template<scalar a, scalar b, IndexAssignment Ix>
    struct simplify_t<ProdExpr<Constant<a>, Constant<b>, Ix>>
    {
        using expr_t = Constant<a * b>;
    };

    template<scalar a, scalar b, IndexAssignment Ix>
    struct simplify_t<SumExpr<Constant<a>, Constant<b>, Ix>>
    {
        using expr_t = Constant<a + b>;
    };

    template<scalar a, NonConstant RLHS, NonConstant LRHS, NonConstant RRHS>
    struct simplify_t<ProdExpr<ProdExpr<Constant<a>, RLHS>, ProdExpr<LRHS, RRHS>>>
    {
        using expr_t = ProdExpr<Constant<a>, ProdExpr<RLHS, ProdExpr<LRHS, RRHS>>>;
    };

    template<scalar a, NonConstant RHS>
    struct simplify_t<FuncExpr<inv_t, ProdExpr<Constant<a>, RHS>>>
    {
        using expr_t = ProdExpr<Constant<a>, FuncExpr<inv_t, RHS>>;
    };

    template<NonConstant E>
    struct simplify_t<ProdExpr<E, FuncExpr<inv_t, E>>>
    {
        using expr_t = OneExpr;
    };

    template<NonConstant E, NonConstant RHS>
    struct simplify_t<ProdExpr<E, FuncExpr<inv_t, ProdExpr<E, RHS>>>>
    {
        using expr_t = FuncExpr<inv_t, RHS>;
    };

    template<NonPower E, NonPower RHS>
    struct simplify_t<ProdExpr<E, ProdExpr<E, RHS>>>
    {
        using expr_t = ProdExpr<PowerExpr<E, 2>, RHS>;
    };

    template<Expr E, int Power>
    struct simplify_t<ProdExpr<E, PowerExpr<E, Power>>>
    {
        using expr_t = PowerExpr<E, Power + 1>;
    };

    template<Expr E, int Power>
    struct simplify_t<ProdExpr<PowerExpr<E, Power>, E>>
    {
        using expr_t = PowerExpr<E, Power + 1>;
    };

    template<NonConstant E, int Power, NonConstant RHS>
    struct simplify_t<ProdExpr<PowerExpr<E, Power>, ProdExpr<E, RHS>>>
    {
        using expr_t = ProdExpr<PowerExpr<E, Power + 1>, RHS>;
    };

    template<NonConstant E, int Power, NonConstant RHS>
    struct simplify_t<ProdExpr<E, ProdExpr<PowerExpr<E, Power>, RHS>>>
    {
        using expr_t = ProdExpr<PowerExpr<E, Power + 1>, RHS>;
    };

    template<NonConstant E, int Power1, int Power2>
    struct simplify_t<ProdExpr<PowerExpr<E, Power1>, PowerExpr<E, Power2>>>
    {
        using expr_t = PowerExpr<E, Power1 + Power2>;
    };

    template<NonConstant E, int Power1, int Power2>
    struct simplify_t<PowerExpr<PowerExpr<E, Power1>, Power2>>
    {
        using expr_t = PowerExpr<E, Power1 * Power2>;
    };

    template<Expr E, int Power>
    struct simplify_t<FuncExpr<inv_t, PowerExpr<E, Power>>>
    {
        using expr_t = PowerExpr<E, -Power>;
    };

    template<Expr E>
    struct simplify_t<FuncExpr<inv_t, E>>
    {
        using expr_t = PowerExpr<E, -1>;
    };

    template<NonConstant LHS, NonConstant RHS>
    struct simplify_t<FuncExpr<inv_t, ProdExpr<LHS, RHS>>>
    {
        using expr_t = ProdExpr<PowerExpr<LHS, -1>, FuncExpr<inv_t, RHS>>;
    };

    template<NonConstant E, int Power1, int Power2, NonConstant RHS>
    struct simplify_t<ProdExpr<PowerExpr<E, Power1>, ProdExpr<PowerExpr<E, Power2>, RHS>>>
    {
        using expr_t = ProdExpr<PowerExpr<E, Power1 + Power2>, RHS>;
    };

    template<scalar a, Expr C, Expr D>
    struct simplify_t<ProdExpr<Constant<a>, SumExpr<C, D>>>
    {
        using expr_t = SumExpr<ProdExpr<Constant<a>, C>, ProdExpr<Constant<a>, D>>;
    };

    template<Expr A, Expr B, Expr C, Expr D>
    struct simplify_t<ProdExpr<SumExpr<A, B>, SumExpr<C, D>>>
    {
        using expr_t = SumExpr<ProdExpr<A, C>, SumExpr<ProdExpr<A, D>, SumExpr<ProdExpr<B, C>, ProdExpr<B, D>>>>;
    };

    template<typename f, Expr ArgExpr, IndexAssignment Ix>
    struct simplify_t<FuncExpr<f, ArgExpr, Ix>>
    {
        using expr_t = FuncExpr<f, typename simplify_t<ArgExpr>::expr_t, Ix>;
    };

    template<Expr Arg>
    struct simplify_t<SumExpr<ProdExpr<FuncExpr<sin_t, Arg>, FuncExpr<sin_t, Arg>>, ProdExpr<FuncExpr<cos_t, Arg>, FuncExpr<cos_t, Arg>>>>
    {
        using expr_t = OneExpr;
    };

    template<Expr Arg>
    struct simplify_t<SumExpr<ProdExpr<FuncExpr<cos_t, Arg>, FuncExpr<cos_t, Arg>>, ProdExpr<FuncExpr<sin_t, Arg>, FuncExpr<sin_t, Arg>>>>
    {
        using expr_t = OneExpr;
    };

    template<Expr Arg, NonConstant S>
    struct simplify_t<SumExpr<ProdExpr<S, ProdExpr<FuncExpr<sin_t, Arg>, FuncExpr<sin_t, Arg>>>, ProdExpr<S, ProdExpr<FuncExpr<cos_t, Arg>, FuncExpr<cos_t, Arg>>>>>
    {
        using expr_t = S;
    };

    template<Expr Arg, NonConstant S, NonConstant S2>
    struct simplify_t<SumExpr<ProdExpr<S2, ProdExpr<S, ProdExpr<FuncExpr<sin_t, Arg>, FuncExpr<sin_t, Arg>>>>, ProdExpr<S2, ProdExpr<S, ProdExpr<FuncExpr<cos_t, Arg>, FuncExpr<cos_t, Arg>>>>>>
    {
        using expr_t = ProdExpr<S2, S>;
    };

    template<VarExp S>
    struct simplify_t<ProdExpr<S, S>>
    {
        using expr_t = PowerExpr<S, 2>;
    };

    template<VarExp S, int Power>
    struct simplify_t<ProdExpr<PowerExpr<S, -Power>, PowerExpr<S, Power>>>
    {
        using expr_t = OneExpr;
    };

    template<int Power, Expr LHS, Expr RHS>
    struct simplify_t<PowerExpr<ProdExpr<LHS, RHS>, Power>>
    {
        using expr_t = ProdExpr<PowerExpr<LHS, Power>, PowerExpr<RHS, Power>>;
    };
 
    template<scalar a>
    struct simplify_t<PowerExpr<Constant<a>, -1>>
    {
        using expr_t = Constant<1.0 / a>;
    };

#if 0
    template<NonConstant A, NonConstant B, NonConstant C, NonConstant D>
    struct simplify_t<SumExpr<SumExpr<A, B>, SumExpr<C, D>>>
    {
        using expr_t = SumExpr<A, SumExpr<B, SumExpr<C, D>>>;
    };
#endif

    template<Expr ExprT>
    constexpr auto simplify_index(ExprT expr)
    {
        if constexpr (has_indexed_expr<ExprT>())
            return expr;
        else if constexpr (std::true_type{})
            return typename unassign_index_t<ExprT>::expr_t{};
    }

    template<Expr ExprT>
    constexpr auto do_simplify(ExprT expr)
    {
        using e1 = decltype(simplify_index(expr));
        using expr_t = simplify_t<e1>::expr_t;

        if constexpr (is_array<expr_t>{})
            return Array<decltype(simplify(typename expr_t::act_t{})), decltype(simplify(typename expr_t::next_t{})), typename expr_t::ix_assign> {};
        else if constexpr (is_sum_expr<expr_t>{})
        {
            if constexpr (is_constant<typename expr_t::lhs_t>{} && is_constant<typename expr_t::rhs_t>{})
                return Constant<expr_t::lhs_t::value + expr_t::rhs_t::value>{};
            else if constexpr (std::is_same_v<typename expr_t::lhs_t, ZeroExpr>)
                return simplify(typename expr_t::rhs_t{});
            else if constexpr (std::is_same_v<typename expr_t::rhs_t, ZeroExpr>)
                return simplify(typename expr_t::lhs_t{});
            else if constexpr (std::is_same_v<typename expr_t::lhs_t, typename expr_t::rhs_t>)
                return ProdExpr<decltype(c<2.0>), decltype(simplify(typename expr_t::lhs_t{})), typename expr_t::ix_assign> {};
            else if constexpr (hash<typename expr_t::lhs_t> > hash<typename expr_t::rhs_t>)
                return simplify(SumExpr<typename expr_t::rhs_t, typename expr_t::lhs_t, typename ExprT::ix_assign>{});
            else if constexpr (std::true_type{})
                return SumExpr<decltype(simplify(typename expr_t::lhs_t{})), decltype(simplify(typename expr_t::rhs_t{})), typename expr_t::ix_assign > {};
        }
        else if constexpr (is_prod_expr<expr_t>{})
        {
            if constexpr (is_constant<typename expr_t::lhs_t>{}&& is_constant<typename expr_t::rhs_t>{})
                return Constant<expr_t::lhs_t::value * expr_t::rhs_t::value>{};
            else if constexpr (std::is_same_v<typename expr_t::lhs_t, ZeroExpr>)
                return ZeroExpr{};
            else if constexpr (std::is_same_v<typename expr_t::rhs_t, ZeroExpr>)
                return ZeroExpr{};
            else if constexpr (std::is_same_v<typename expr_t::lhs_t, OneExpr>)
                return simplify(typename expr_t::rhs_t{});
            else if constexpr (std::is_same_v<typename expr_t::rhs_t, OneExpr>)
                return simplify(typename expr_t::lhs_t{});
            else if constexpr (hash<typename expr_t::lhs_t> > hash<typename expr_t::rhs_t>)
                return simplify(ProdExpr<typename expr_t::rhs_t, typename expr_t::lhs_t, typename ExprT::ix_assign>{});
            else if constexpr (is_prod_expr<typename expr_t::rhs_t>{})
            {
                if constexpr (hash<typename expr_t::lhs_t> > hash<typename expr_t::rhs_t::lhs_t>)
                    return ProdExpr<typename expr_t::rhs_t::lhs_t, ProdExpr<typename expr_t::lhs_t, typename expr_t::rhs_t::rhs_t, typename expr_t::ix_assign>, typename expr_t::ix_assign>{};
                else if constexpr (std::true_type{})
                    return ProdExpr<decltype(simplify(typename expr_t::lhs_t{})), decltype(simplify(typename expr_t::rhs_t{})), typename expr_t::ix_assign > {};
            }
            else if constexpr (std::true_type{})
                return ProdExpr<decltype(simplify(typename expr_t::lhs_t{})), decltype(simplify(typename expr_t::rhs_t{})), typename expr_t::ix_assign > {};
        }
        else if constexpr (is_func_expr<expr_t>{})
        {
            return FuncExpr<typename expr_t::f_t, decltype(simplify(typename expr_t::arg_t{})), typename expr_t::ix_assign > {};
        }
        else if constexpr (is_power_expr<expr_t>{})
        {
            return PowerExpr<decltype(simplify(typename expr_t::expr_t{})), expr_t::power, typename expr_t::ix_assign > {};
        }
        else if constexpr (std::true_type{})
            return expr_t{};
    }

    template<Expr ExprT>
    auto simplify(ExprT expr)
    {
        using R = decltype(do_simplify(expr));
        if constexpr (std::is_same_v<ExprT, R>)
            return expr;
        else if constexpr (std::true_type{})
            return simplify(R{});
    }

    template<Expr ExprT>
    constexpr bool all_index_assigned()
    {
        if constexpr (is_array<ExprT>{})
            return all_index_assigned<typename ExprT::act_t>() && all_index_assigned<typename ExprT::next_t>();
        else if constexpr (is_indexed_expr<ExprT>{})
            return all_index_assigned<typename ExprT::expr_t>();
        else if constexpr (is_sum_expr<ExprT>{} || is_prod_expr<ExprT>{})
            return all_index_assigned<typename ExprT::lhs_t>() && all_index_assigned<typename ExprT::rhs_t>();
        else if constexpr (is_derivative_expr<ExprT>{})
            return all_index_assigned<typename ExprT::f_t>() && all_index_assigned<typename ExprT::v_t>();
        else if constexpr (is_func_expr<ExprT>{})
            return all_index_assigned<typename ExprT::f_t>() && all_index_assigned<typename ExprT::arg_t>();
        else if constexpr (is_var_expr<ExprT>{})
            return is_index_assigned<typename ExprT::ix_assign, ExprT::V.ix>();
        else if (std::true_type{})
        {
            debug<ExprT>{};
            return false;
        }
    }

    template<Expr ExprT>
    constexpr bool has_indexed_expr()
    {
        if constexpr (is_indexed_expr<ExprT>{})
            return true;
        else if constexpr (std::is_same_v<ExprT, nothing>)
            return false;
        else if constexpr (is_var_expr<ExprT>{})
            return false;
        else if constexpr (is_constant<ExprT>{})
            return false;
        else if constexpr (is_array<ExprT>{})
            return has_indexed_expr<typename ExprT::act_t>() || has_indexed_expr<typename ExprT::next_t>();
        else if constexpr (is_sum_expr<ExprT>{} || is_prod_expr<ExprT>{})
            return has_indexed_expr<typename ExprT::lhs_t>() || has_indexed_expr<typename ExprT::rhs_t>();
        else if constexpr (is_derivative_expr<ExprT>{})
            return has_indexed_expr<typename ExprT::f_t>() || has_indexed_expr<typename ExprT::v_t>();
        else if constexpr (is_func_expr<ExprT>{})
            return has_indexed_expr<typename ExprT::arg_t>();
        else if constexpr (is_power_expr<ExprT>{})
            return has_indexed_expr<typename ExprT::expr_t>();
        else if constexpr (std::true_type{})
        {
            debug<ExprT>{};
            return true;
        }
    }
}
