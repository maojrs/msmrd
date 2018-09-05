//Based on ReaDDy vec3 class

#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <sstream>
#include <ostream>

//namespace py = pybind11;

namespace detail {
    template<typename T>
    using is_arithmetic_type = std::enable_if_t<std::is_arithmetic<T>::value, int>;
}

template<typename scalar=double>
class vec3 {
public:
    using data_arr = std::array<scalar, 3>;

    union {
        struct {
            scalar x, y, z;
        };
        data_arr data;
    };

    vec3(): vec3(0, 0, 0){};

    vec3(scalar x, scalar y, scalar z): x(x), y(y), z(z) {};

    vec3(std::vector<double> &v): x(v[0]), y(v[1]), z(v[2]) {};

    //explicit vec3(const data_arr &xyz) : data(xyz) {};
    vec3 &operator+=(const vec3 &rhs) {
        std::transform(data.begin(), data.end(), rhs.data.begin(), data.begin(), std::plus<scalar>());
        return *this;
    };

    vec3 &operator-=(const vec3 &rhs) {
        std::transform(data.begin(), data.end(), rhs.data.begin(), data.begin(), std::minus<scalar>());
        return *this;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    vec3 &operator+=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind1st(std::plus<scalar>(), a));
        return *this;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    vec3 &operator-=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind2nd(std::minus<scalar>(), a));
        return *this;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    vec3 &operator*=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind1st(std::multiplies<scalar>(), a));
        return *this;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    vec3 &operator/=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind2nd(std::divides<scalar>(), a));
        return *this;
    };

    vec3 cross(const vec3 &other) const {
        return {
                data[1] * other.data[2] - data[2] * other.data[1],
                data[2] * other.data[0] - data[0] * other.data[2],
                data[0] * other.data[1] - data[1] * other.data[0]
        };
    };

    scalar norm() const {
        return std::sqrt(normSquared());
    };

    scalar normSquared() const {
        return std::inner_product(data.begin(), data.end(), data.begin(), static_cast<scalar>(0));
    };

    scalar infnorm() const {
        scalar inorm = 0;
        for (auto i = 0; i < 3; ++i) {
            inorm += std::abs(data[i]);
        }
        return inorm;
    };

    scalar operator[](std::size_t i) const {
        return data.at(i);
    };

    scalar &operator[](std::size_t i) {
        return data.at(i);
    };

    vec3 &invertElementWise() {
        for (auto i = 0; i < 3; ++i) {
            data[i] = static_cast<scalar>(1.) / data[i];
        }
        return *this;
    };

    bool operator==(const vec3 &rhs) const {
        return data == rhs.data;
    };

    bool operator!=(const vec3 &rhs) const {
        return data != rhs.data;
    };

//    bool almostEquals(const vec3 &rhs) const {
//        bool result {true};
//        for(std::uint8_t i = 0; i < 3; ++i) {
//            const auto fp1 = fp::FloatingPoint<float>(data[i]);
//            const auto fp2 = fp::FloatingPoint<float>(rhs.data[i]);
//            auto x = fp::FloatingPoint<float>::DistanceBetweenSignAndMagnitudeNumbers(fp1.bits(), fp2.bits());
//            result &= fp1.AlmostEquals(fp2);
//        }
//        return result;
//    }

    friend std::ostream &operator<<(std::ostream &os, const vec3 &vec) {
        os << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
        return os;
    };

    friend vec3 operator+(vec3 lhs, const vec3 &rhs) {
        lhs += rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend vec3 operator+(vec3 lhs, arithmetic rhs) {
        lhs += rhs;
        return lhs;
    };

    friend vec3 operator-(vec3 lhs, const vec3 &rhs) {
        lhs -= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend vec3 operator-(vec3 lhs, arithmetic rhs) {
        lhs -= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend vec3 operator/(vec3 lhs, arithmetic rhs) {
        lhs /= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend vec3 operator*(vec3 lhs, const arithmetic rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend bool operator>=(const vec3 &lhs, const vec3 &rhs) {
        bool result{true};
        for (auto i = 0; i < 3; ++i) { result &= lhs[i] >= rhs[i]; }
        return result;
    };

    friend bool operator<=(const vec3 &lhs, const vec3 &rhs) {
        bool result{true};
        for (auto i = 0; i < 3; ++i) { result &= lhs[i] <= rhs[i]; }
        return result;
    };

    friend bool operator>(const vec3 &lhs, const vec3 &rhs) {
        return !(lhs <= rhs);
    };

    friend bool operator<(const vec3 &lhs, const vec3 &rhs) {
        return !(lhs >= rhs);
    };

    friend scalar operator*(const vec3 &lhs, const vec3 &rhs) {
        return std::inner_product(lhs.data.begin(), lhs.data.end(), rhs.data.begin(), static_cast<scalar>(0));
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend vec3 operator*(const arithmetic rhs, vec3 lhs) {
        lhs *= rhs;
        return lhs;
    }

};