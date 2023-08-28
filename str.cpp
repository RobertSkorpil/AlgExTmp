#include "str.h"

std::wstring_convert<codecvt_u16, char16_t> convert;

std::string to_str(const std::u16string& s)
{
    return convert.to_bytes(s);
}

std::u16string to_u16str(const std::string& s)
{
    return convert.from_bytes(s);
}
