#include <type_traits>
#include <string>
#include <optional>
#include <algorithm>
#include <math.h>
#include "str.h"

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
        index ix{};

        constexpr var(const var&) = default;
        constexpr var(char16_t n) : n{ n } {}
        constexpr var(char16_t n, char16_t ix) : n{ n }, ix{ ix } {}

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
    concept NonConstant = requires { T::non_const; }&& Expr<T>;

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
    
    template<index Ix, size_t Val, IndexAssignment Next>
    struct index_assignment
    {
        static constexpr index ix { Ix };
        static constexpr size_t val { Val };
        using next = Next;
    };

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

    template<Expr ExprT, index ix, IndexAssignment Ix>
    struct IndexedExpr;

    using func_single_t = scalar(scalar);

    template<scalar c, IndexAssignment Ix = no_index>
    struct Constant
    {
        using ix_assign = Ix;
        template<typename Valuation>
        static auto eval(const Valuation& ctx)
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

    template<var v, IndexAssignment Ix = no_index>
    struct VarExpr
    {
        static constexpr bool is_var{ true };
        static constexpr bool non_const{ true };
        static constexpr var V { v };
        using ix_assign = Ix;

        template<typename Valuation>
        static auto eval(const Valuation& ctx)
        {
            return ctx.var_value(v, ix_val<Ix>(v.ix));
        }

        template<VarExp D>
        constexpr static auto diff()
        {
            constexpr auto d{ D::V };
            if constexpr (d == v)
                return OneExpr{};
            else if constexpr (d.n == v.n)
            {
                constexpr auto ix1{ ix_val<Ix>(v.ix) };
                constexpr auto ix2{ ix_val<typename D::ix_assign>(D::V.ix) };
                if constexpr (ix1 && ix2 && *ix1 == *ix2)
                    return OneExpr{};
                else
                    return ZeroExpr{};
            }
            else if constexpr (d != v)
                return ZeroExpr{};
            else
                return ZeroExpr{};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            auto ix{ ix_val<Ix>(v.ix) };
            return std::u16string{ v.n } + (v.ix.n ? (ix ? std::u16string{ wchar_t(u'₀' + *ix) } : std::u16string{ v.ix.n }) : u"");
        }

        template<Expr RHS>
        constexpr auto operator =(const RHS& rhs);
    };

    template<typename T>
    struct is_var_expr : std::false_type {};

    template<var v, IndexAssignment Ix>
    struct is_var_expr<VarExpr<v, Ix>> : std::true_type{};

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    struct Assignment
    {
        static constexpr bool is_action{ true };
        using ix_assign = Ix;

        template<typename Valuation>
        static void eval(const Valuation& ctx)
        {
            ctx.set_var_value(LHS::V, ix_val<LHS::Ix>(LHS::V.ix), RHS::eval(ctx));
        }

        static std::u16string to_str()
        {
            return std::remove_reference_t<LHS>::to_str() + u" = " + RHS::to_str();
        }
    };

    struct no_action
    {
        using act_t = no_action;
        using next_t = no_action;

        static constexpr bool is_action{ true };

        template<typename Valuation>
        static void eval(const Valuation& ctx)
        {}


        static std::u16string to_str()
        {
            return {};
        }
    };

    template<Action Act, Action Next = no_action, IndexAssignment Ix = no_index>
    struct Array
    {
        static constexpr bool is_action{ true };
        using act_t = Act;
        using next_t = Next;
        using ix_assign = Ix;

        template<typename Valuation>
        static void eval(const Valuation& ctx)
        {
            Act::eval(ctx);
            Next::eval(ctx);
        }

        static std::u16string to_str()
        {
            return Act::to_str() + u"\n" + Next::to_str();
        }
    };

    template<typename T>
    struct is_array : std::false_type {};

    template<Action Act, Action Next, IndexAssignment Ix>
    struct is_array<Array<Act, Next, Ix>> : std::true_type{};

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
        static constexpr int non_const{};
        using lhs_t = LHS;
        using rhs_t = RHS;
        using ix_assign = Ix;

        template<typename Valuation>
        static auto eval(const Valuation& ctx)
        {
            return LHS::eval(ctx) + RHS::eval(ctx);
        }

        template<VarExp D>
        static constexpr SumExpr<typename LHS::template diff_t<D>, typename RHS::template diff_t<D>, Ix> diff()
        {
            return {};
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
        static auto eval(const Valuation& ctx)
        {
            return LHS::eval(ctx) * RHS::eval(ctx);
        }

        template<VarExp D>
        static constexpr SumExpr<ProdExpr<LHS, typename RHS::template diff_t<D>>, ProdExpr<typename LHS::template diff_t<D>, RHS>> diff()
        {
            return {};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            return u"(" + LHS::to_str() + u" * " + RHS::to_str() + u")";
        }
    };

    template<typename T>
    struct is_prod_expr : std::false_type {};

    template<Expr LHS, Expr RHS, IndexAssignment Ix>
    struct is_prod_expr<ProdExpr<LHS, RHS, Ix>> : std::true_type {};

    template<typename f, Expr ArgExpr, IndexAssignment Ix = no_index>
    struct FuncExpr;

    template<func_single_t f, StringLiteral name, IndexAssignment Ix = no_index>
    struct Func
    {
        static double eval(double arg)
        {
            return f(arg);
        }

        template<Expr ArgExpr>
        constexpr auto operator()(const ArgExpr& arg)
        {
            return FuncExpr<Func, ArgExpr, Ix>{};
        }

        static std::u16string to_str()
        {
            return name.to_str();
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

    template<Expr F, VarExp V, IndexAssignment Ix = no_index>
    struct DerivativeExpr
    {
        static constexpr int non_const{};
        using f_t = F;
        using v_t = V;
        using ix_assign = Ix;

        template<typename Valuation>
        static auto eval(const Valuation& ctx)
        {
            return std::numeric_limits<scalar>::quiet_NaN();
        }

        constexpr static auto derive()
        {
            if constexpr (!V::V.ix)
                return typename F::template diff_t<V>{};
            else if constexpr (all_index_assigned<F>() && is_index_assigned<typename V::ix_assign, V::V.ix>())
                return typename F::template diff_t<V>{};
            else if (std::true_type{})
                return DerivativeExpr{};
        }

        static std::u16string to_str()
        {
            return u" d(" + F::to_str() + u")/d(" + V::to_str() + u") ";
        }
    };

    template<typename T>
    struct is_derivative_expr : std::false_type {};

    template<Expr F, VarExp V, IndexAssignment Ix>
    struct is_derivative_expr<DerivativeExpr<F, V, Ix>> : std::true_type {};

    template<DifferentialLike F, DifferentialLike V>
    constexpr auto operator /(const F&, const V&)
    {
        using DE = DerivativeExpr<typename F::expr_t, typename V::expr_t>;
        return DE::derive();
    }

    template<typename f, Expr ArgExpr, IndexAssignment Ix>
    struct func_derivative {};

    using sin_t = Func<sin, "sin">; sin_t Sin{};
    using cos_t = Func<cos, "cos">; cos_t Cos{};

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

    template<typename F, Expr ArgExpr, IndexAssignment Ix>
    struct FuncExpr
    {
        static constexpr int non_const{};
        using f_t = F;
        using arg_t = ArgExpr;
        using ix_assign = Ix;

        template<typename Valuation>
        static auto eval(const Valuation& ctx)
        {
            return F::eval(ArgExpr::eval(ctx));
        }

        template<VarExp D>
        static constexpr
            ProdExpr< typename func_derivative<F, ArgExpr, Ix>::expr_t, typename ArgExpr::template diff_t<D>, Ix> diff()
        {
            return {};
        }

        template<VarExp D>
        using diff_t = decltype(diff<D>());

        static std::u16string to_str()
        {
            return F::to_str() + u"(" + ArgExpr::to_str() + u")";
        }
    };

    template<Expr ExprT, index ix, IndexAssignment Ix = no_index>
    struct IndexedExpr
    {
        using expr_t = ExprT;        
        using ix_assign = Ix;
        static constexpr index expr_ix = ix;

        template<typename Valuation>
        static auto eval(const Valuation& ctx)
        {
            return std::numeric_limits<scalar>::quiet_NaN();
        }

        template<Expr RHS>
        constexpr auto operator =(const RHS& rhs)
        {
            return Assignment<decltype(*this), RHS, Ix>{};
        }

        static std::u16string to_str()
        {
            return ExprT::to_str() + u"[" + ix.n + u"]";
        }
    };

    template<typename T>
    struct is_indexed_expr : std::false_type {};

    template<Expr ExprT, index ix, IndexAssignment Ix>
    struct is_indexed_expr<IndexedExpr<ExprT, ix, Ix>> : std::true_type {};

    template<index ix, Expr ExprT>
    constexpr auto Ix(ExprT expr)
    {
        return IndexedExpr<ExprT, ix>{};
    }

    template<typename T>
    struct is_func_expr : std::false_type {};

    template<typename F, Expr ArgExpr, IndexAssignment Ix>
    struct is_func_expr<FuncExpr<F, ArgExpr, Ix>> : std::true_type {};

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    auto operator +(LHS lhs, RHS rhs)
    {
        return SumExpr<LHS, RHS, Ix>{};
    }

    template<Expr LHS, Expr RHS, IndexAssignment Ix = no_index>
    auto operator *(LHS lhs, RHS rhs)
    {
        return ProdExpr<LHS, RHS, Ix>{};
    }

    template<scalar val>
    constexpr Constant<val> c{};

    template<index ix, size_t val, Expr ExprT>
    struct assign_index_t
    {
    };

    template<index ix, size_t val, IndexAssignment Next, Expr ExprT>
    struct assign_index_t<ix, val, IndexedExpr<ExprT, ix, Next>>
    {
        using expr_t = IndexedExpr<ExprT, ix, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr ExprT, index expr_ix>
    struct assign_index_t<ix, val, IndexedExpr<ExprT, expr_ix, Next>>
    {
        using expr_t = IndexedExpr<typename assign_index_t<ix, val, ExprT>::expr_t, expr_ix, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, scalar c>
    struct assign_index_t<ix, val, Constant<c, Next>>
    {
        using expr_t = Constant<c, Next>;
    };

    template<index ix, size_t val>
    struct assign_index_t<ix, val, no_action>
    {
        using expr_t = no_action;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, SumExpr<LHS, RHS, Next>>;

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, ProdExpr<LHS, RHS, Next>>;

    template<index ix, size_t val, IndexAssignment Next, typename F, Expr Arg>
    struct assign_index_t<ix, val, FuncExpr<F, Arg, Next>>;

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, DerivativeExpr<LHS, RHS, Next>>;

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, Assignment<LHS&, RHS, Next>>;

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, Assignment<LHS, RHS, Next>>;

    template<index ix, size_t val, IndexAssignment Next, var v>
    struct assign_index_t<ix, val, VarExpr<v, Next>>;

    template<index ix, size_t val, IndexAssignment Next, Action Act, Action NextAct>
    struct assign_index_t<ix, val, Array<Act, NextAct, Next>>;

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, SumExpr<LHS, RHS, Next>>
    {
        using expr_t = SumExpr<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, ProdExpr<LHS, RHS, Next>>
    {
        using expr_t = ProdExpr<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, DerivativeExpr<LHS, RHS, Next>>
    {
        using expr_t = DerivativeExpr<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, Assignment<LHS&, RHS, Next>>
    {
        using expr_t = Assignment<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, Expr LHS, Expr RHS>
    struct assign_index_t<ix, val, Assignment<LHS, RHS, Next>>
    {
        using expr_t = Assignment<typename assign_index_t<ix, val, LHS>::expr_t, typename assign_index_t<ix, val, RHS>::expr_t, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, var v>
    struct assign_index_t<ix, val, VarExpr<v, Next>>
    {
        using expr_t = VarExpr<v, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, Action Act, Action NextAct>
    struct assign_index_t<ix, val, Array<Act, NextAct, Next>>
    {
        using expr_t = Array<typename assign_index_t<ix, val, Act>::expr_t, typename assign_index_t<ix, val, NextAct>::expr_t, index_assignment<ix, val, Next>>;
    };

    template<index ix, size_t val, IndexAssignment Next, typename F, Expr Arg>
    struct assign_index_t<ix, val, FuncExpr<F, Arg, Next>>
    {
        using expr_t = FuncExpr<F, typename assign_index_t<ix, val, Arg>::expr_t, index_assignment<ix, val, Next>>;
    };

    template<Expr ExprT>
    constexpr auto expand(ExprT expr);

    template<index ix, size_t val, size_t range, Expr ExprT>
    auto sum_sumand(ExprT expr);

    template<index ix, size_t val, size_t range, Expr ExprT>
    auto sum_sumand(ExprT expr)
    {
        if constexpr (val == range)
            return ZeroExpr{};
        else if constexpr (val == range - 1)
            return typename assign_index_t<ix, val, decltype(expr)>::expr_t{};
        else
            return SumExpr<typename assign_index_t<ix, val, decltype(expr)>::expr_t, decltype(sum_sumand<ix, val + 1, range, ExprT>(expr)), index_assignment<ix, val, typename ExprT::ix_assign>>{};
    }

    template<index ix, size_t range, Expr ExprT>
    auto Sum(ExprT expr)
    {
        return simplify(expand(sum_sumand<ix, 0, range>(simplify(expand(expr)))));
    }

    template<index ix, size_t val, size_t range, Expr ExprT>
    constexpr auto make_for(ExprT expr);

    template<index ix, size_t val, size_t range, Expr ExprT>
    constexpr auto make_for(ExprT expr)
    {
        if constexpr (val == range)
            return no_action{};
        else
            return Array<typename assign_index_t<ix, val, decltype(expr)>::expr_t, decltype(make_for<ix, val + 1, range, ExprT>(expr))> {};
    }

    template<index ix, size_t range, Expr ExprT>
    constexpr auto For(ExprT expr)
    {
        return simplify(expand(make_for<ix, 0, range>(simplify(expand(expr)))));
    }

    template<size_t i, size_t v, Expr Main, Expr ExprT>
    constexpr auto do_index_seq(ExprT expr);

    template<typename T>
    struct debug;

    template<size_t v>
    struct debug_size_t;

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

    template<Expr ExprT>
    constexpr auto do_index(ExprT expr)
    {
        if constexpr (is_array<typename ExprT::expr_t>{})
        {
            constexpr auto v{ ix_val_static<typename ExprT::ix_assign, ExprT::expr_ix>() };
            return do_index_seq<size_t(0), v, ExprT>(typename ExprT::expr_t{});
        }
        else if constexpr (std::true_type{})
            return expr;
    }

    template<Expr ExprT>
    constexpr auto do_expand(ExprT expr);

    template<Expr ExprT>
    constexpr auto do_expand(ExprT expr)
    {
        if constexpr (is_array<ExprT>{})
            return Array<decltype(do_expand(typename ExprT::act_t{})), decltype(do_expand(typename ExprT::next_t{})), typename ExprT::ix_assign > {};
        else if constexpr (is_indexed_expr<ExprT>{})
        {
            if constexpr (is_index_assigned<typename ExprT::ix_assign, ExprT::expr_ix>())
                return do_index(expr);
            else if constexpr (std::true_type{})
                return expr;
        }
        else if constexpr (is_sum_expr<ExprT>{})
            return SumExpr<decltype(do_expand(typename ExprT::lhs_t{})), decltype(do_expand(typename ExprT::rhs_t{})), typename ExprT::ix_assign > {};
        else if constexpr (is_prod_expr<ExprT>{})
            return ProdExpr<decltype(do_expand(typename ExprT::lhs_t{})), decltype(do_expand(typename ExprT::rhs_t{})), typename ExprT::ix_assign > {};
        else if constexpr (is_derivative_expr<ExprT>{})
            return ExprT::derive();
        else if constexpr (is_func_expr<ExprT>{})
            return FuncExpr<decltype(do_expand(typename ExprT::f_t{})), decltype(do_expand(typename ExprT::arg_t{})), typename ExprT::ix_assign > {};
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

    template<Expr ProdLHS, Expr ProdRHS, IndexAssignment Ix>
    struct simplify_t<SumExpr<ProdExpr<ProdLHS, ProdRHS, Ix>, ProdExpr<ProdRHS, ProdLHS, Ix>, Ix>>
    {
        using slhs = simplify_t<ProdLHS>::expr_t;
        using srhs = simplify_t<ProdRHS>::expr_t;
        using expr_t = simplify_t<SumExpr<ProdExpr<slhs, srhs, Ix>, ProdExpr<slhs, srhs, Ix>, Ix>>::expr_t;
    };

    template<Expr A, Expr B, IndexAssignment Ix>
    struct simplify_t<ProdExpr<SumExpr<A, B, Ix>, SumExpr<A, B, Ix>, Ix>>
    {
        using sa = simplify_t<A>::expr_t;
        using sb = simplify_t<B>::expr_t;

        using expr_t = SumExpr<SumExpr<ProdExpr<sa, sa, Ix>, ProdExpr<sb, sb, Ix>>, ProdExpr<Constant<2.0>, SumExpr<sa, sb, Ix>, Ix>, Ix>;
    };

    template<Expr A, Expr B, Expr C, Expr D, IndexAssignment Ix>
    struct simplify_t<ProdExpr<SumExpr<A, B, Ix>, SumExpr<C, D, Ix>, Ix>>
    {
        using sa = simplify_t<A>::expr_t;
        using sb = simplify_t<B>::expr_t;
        using sc = simplify_t<C>::expr_t;
        using sd = simplify_t<D>::expr_t;

        using expr_t = SumExpr<SumExpr<ProdExpr<sa, sc, Ix>, ProdExpr<sa, sd, Ix>>, SumExpr<ProdExpr<sb, sc, Ix>, ProdExpr<sb, sd, Ix>, Ix>, Ix>;
    };

    template<typename f, Expr ArgExpr, IndexAssignment Ix>
    struct simplify_t<FuncExpr<f, ArgExpr, Ix>>
    {
        using expr_t = FuncExpr<f, typename simplify_t<ArgExpr>::expr_t, Ix>;
    };

    template<Expr ExprT>
    constexpr auto do_simplify(ExprT expr)
    {
        using expr_t = ExprT;// simplify_t<ExprT>;

        if constexpr (is_array<expr_t>{})
            return Array<decltype(do_simplify(typename expr_t::act_t{})), decltype(do_simplify(typename expr_t::next_t{})), typename expr_t::ix_assign> {};
        else if constexpr (is_sum_expr<expr_t>{})
        {
            if constexpr (std::is_same_v<typename expr_t::lhs_t, ZeroExpr>)
                return do_simplify(typename expr_t::rhs_t{});
            else if constexpr (std::is_same_v<typename expr_t::rhs_t, ZeroExpr>)
                return do_simplify(typename expr_t::lhs_t{});
            else if constexpr (std::is_same_v<typename expr_t::lhs_t, typename expr_t::rhs_t>)
                return ProdExpr<decltype(c<2.0>), decltype(do_simplify(typename expr_t::lhs_t{})), typename expr_t::ix_assign> {};
            else if (std::true_type{})
                return SumExpr<decltype(simplify(typename expr_t::lhs_t{})), decltype(do_simplify(typename expr_t::rhs_t{})), typename expr_t::ix_assign > {};
        }
        else if constexpr (is_prod_expr<expr_t>{})
        {
            if constexpr (std::is_same_v<typename expr_t::lhs_t, ZeroExpr>)
                return ZeroExpr{};
            else if constexpr (std::is_same_v<typename expr_t::rhs_t, ZeroExpr>)
                return ZeroExpr{};
            else if constexpr (std::is_same_v<typename expr_t::lhs_t, OneExpr>)
                return do_simplify(typename expr_t::rhs_t{});
            else if constexpr (std::is_same_v<typename expr_t::rhs_t, OneExpr>)
                return do_simplify(typename expr_t::lhs_t{});
            else if (std::true_type{})
                return ProdExpr<decltype(simplify(typename expr_t::lhs_t{})), decltype(do_simplify(typename expr_t::rhs_t{})), typename expr_t::ix_assign > {};
        }
        else if constexpr (std::true_type{})
            return expr;
    }

    template<Expr ExprT>
    auto simplify(ExprT expr)
    {
        using R = decltype(do_simplify(expr));
        if constexpr (std::is_same_v<ExprT, R>)
            return expr;
        else
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
            return false;
    }

#if 0
    template<var d, Expr ExprT>
    auto diff(ExprT expr)
    {
        return simplify(typename ExprT::template diff_t<d>{});
    }
#endif
}
