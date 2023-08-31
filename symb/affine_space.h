#pragma once
#include <array>
#include <cstdlib>
#include <algorithm>

template<size_t D>
struct affine_space
{
    template<typename T>
    using arr_t = std::array<T, 4>;
    using coords_t = arr_t<double>;
    
    using tensor2 = arr_t<coords_t>;
    using tensor3 = arr_t<arr_t<coords_t>>;

    struct vec;
    struct point
    {
        coords_t coords;

        double* data() {
            return coords.data();
        }

        double operator[](size_t i) const {
            return coords[i];
        }

        double& operator[](size_t i) {
            return coords[i];
        }

        vec operator -(const point& b) const;
        point operator +(const vec& v) const;
    };

    struct vec
    {
        coords_t coords;

        double* data() {
            return coords.data();
        }

        double operator[](size_t i) const {
            return coords[i];
        }

        double& operator[](size_t i) {
            return coords[i];
        }

        vec operator *(double f) const
        {
            vec r;
            std::transform(coords.begin(), coords.end(), r.coords.begin(), [f](double a) { return a * f; });
            return r;
        }

        vec operator +(const vec& v) const
        {
            vec r;
            std::transform(coords.begin(), coords.end(), v.coords.begin(), r.coords.begin(), [](double a, double b) { return a + b; });
            return r;
        }
    };
};

template<size_t D>
inline typename affine_space<D>::vec affine_space<D>::point::operator -(const typename affine_space<D>::point& b) const
{
    typename affine_space<D>::vec r;
    std::transform(coords.begin(), coords.end(), b.begin(), r.begin(), [](double a, double b) { return a - b; });
    return r;
}

template<size_t D>
inline typename affine_space<D>::point affine_space<D>::point::operator +(const typename affine_space<D>::vec& v) const
{
    typename affine_space<D>::point r;
    std::transform(coords.begin(), coords.end(), v.coords.begin(), r.coords.begin(), [](double a, double b) { return a + b; });
    return r;
}
