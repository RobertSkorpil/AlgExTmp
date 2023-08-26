#include <iostream>
#include <locale>
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#include "symb.h"
#include "str.h"

using namespace std;

int main()
{
    SetConsoleOutputCP(CP_UTF8);
    using namespace symb;
    VarExpr<'X'> x;
    VarExpr<'Y'> y;
    VarExpr<'Z'> z;
    VarExpr<var{'Y', 'u'}> yu;
    VarExpr<var{'X', 'u'}> xu;
    VarExpr<var{'X', 'v'}> xv;

    auto dl{ For<'u', 2>(d(xu * xu * xu) / d(xu)) };
    auto dk{ For<'v', 2>(Ix<'v'>(dl)) };

    std::cout << to_str(dl.to_str()) << '\n';

    std::cout << to_str(dk.to_str()) << '\n';
}
