#pragma once

#define _SILENCE_CXX17_CODECVT_HEADER_DEPRECATION_WARNING
#define _SILENCE_CXX20_CODECVT_FACETS_DEPRECATION_WARNING
#include <string>
#include <locale>

struct codecvt_u16 : std::codecvt<char16_t, char, std::mbstate_t>
{
  codecvt_u16(std::size_t refs = 0) : codecvt{refs} {}
};

extern std::wstring_convert<codecvt_u16, char16_t> convert;

std::string to_str(const std::u16string& s);

std::u16string to_u16str(const std::string& s);
