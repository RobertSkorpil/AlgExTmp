#include "str.h"

std::wstring_convert<std::codecvt<char16_t, char, std::mbstate_t>, char16_t> convert;

std::string to_str(const std::u16string& s)
{
    return convert.to_bytes(s);
}

std::u16string to_u16str(const std::string& s)
{
    return convert.from_bytes(s);
}
