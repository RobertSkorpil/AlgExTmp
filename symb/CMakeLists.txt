add_library(symb symb.h symb.cpp "affine_space.h")

set_property(TARGET symb PROPERTY CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
target_include_directories(symb INTERFACE .)