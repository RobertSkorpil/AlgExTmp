#include "symb.h"

namespace symb
{
    sin_t Sin{};
    cos_t Cos{};
    inv_t Inv{};
    sqrt_t Sqrt{};

    scalar do_inverse(scalar x)
    {
        return 1.0 / x;
    }
}