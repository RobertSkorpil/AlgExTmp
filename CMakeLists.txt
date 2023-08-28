# CMakeList.txt : CMake project for AlgExTmp, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.25)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("AlgExTmp")
add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

# Add source to this project's executable.
add_executable (AlgExTmp "AlgExTmp.cpp"  "symb.h" "str.h" "str.cpp" "symb.cpp" "kerr.h" "kerr.cpp" "riemann.h" "matrix.h")

set_property(TARGET AlgExTmp PROPERTY CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

# TODO: Add tests and install targets if needed.
