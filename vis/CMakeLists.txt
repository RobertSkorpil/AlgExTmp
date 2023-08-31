add_library(vis vis.h vis.cpp)

target_include_directories(vis INTERFACE .)

target_link_libraries(vis PRIVATE glfw GLEW::GLEW)