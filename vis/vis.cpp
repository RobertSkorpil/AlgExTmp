#include "vis.h"
#include <gl/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string_view>

using namespace vis;

namespace
{
    void load_shader(uint32_t shaderId, std::string_view sourceFile)
    {
        std::ifstream ifs(std::string{ sourceFile });
        ifs.exceptions(std::ios_base::failbit | std::ios_base::badbit);
        std::string src(std::istreambuf_iterator<char>{ifs}, std::istreambuf_iterator<char>{});

        auto src_data = src.data();
        auto src_len = GLint(src.length());
        glShaderSource(shaderId, 1, &src_data, &src_len);
        glCompileShader(shaderId);

        int success;
        glGetShaderiv(shaderId, GL_COMPILE_STATUS, &success);

        if (!success)
        {
            char info[512];
            glGetShaderInfoLog(shaderId, sizeof(info), nullptr, info);
            std::cerr << info << std::endl;
            throw std::runtime_error("shader syntax error");
        }
    }

    void link_shaders(uint32_t program, std::initializer_list<uint32_t> shaders)
    {
        for (auto id : shaders)
            glAttachShader(program, id);

        glLinkProgram(program);

        int success;
        glGetProgramiv(program, GL_LINK_STATUS, &success);

        if (!success)
        {
            char info[512];
            glGetProgramInfoLog(program, sizeof(info), nullptr, info);
            std::cerr << info << std::endl;
            throw std::runtime_error("shader link error");
        }
    }

    GLFWwindow* window;
}

void vis::init()
{
    glfwInit();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    window = glfwCreateWindow(1200, 1200, "", nullptr, nullptr);
    glfwMakeContextCurrent(window);

    glewInit();

    glEnable(GL_ARB_separate_shader_objects);

    auto vtx_shader = glCreateShader(GL_VERTEX_SHADER);
    load_shader(vtx_shader, "vertex.glsl");

    auto frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    load_shader(frag_shader, "fragment.glsl");

    auto program = glCreateProgram();
    link_shaders(program, { vtx_shader, frag_shader });

    uint32_t vbo;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    uint32_t vao, vao2;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glVertexAttribPointer(0, 4, GL_DOUBLE, GL_FALSE, sizeof(point), 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 4, GL_DOUBLE, GL_FALSE, sizeof(point), (void *)sizeof(vec));
    glEnableVertexAttribArray(1);

    glUseProgram(program);

    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(3);
    glClearColor(0, 0, 0, 0);
}

void vis::swap()
{
    glfwSwapBuffers(window);
    glClear(GL_COLOR_BUFFER_BIT);
    glfwPollEvents();
}

void vis::draw_point(const point &p)
{
    draw_points(&p, 1);
}

void vis::draw_points(const point* p, size_t n)
{
    glBufferData(GL_ARRAY_BUFFER, sizeof(point) * n, p, GL_DYNAMIC_DRAW);
    glDrawArrays(GL_POINTS, 0, n);
}
