#pragma once

#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <type_traits>
#include <sstream>
#include <ostream>
#include "vec3.hpp"

// Quaternion class with operator overloading
template<typename scalar=double>
class quaternion {
public:
    using data_arr = std::array<scalar, 4>;

    union {
        struct {
            scalar a, b, c, d;
        };
        struct {
            scalar re;
            vec3<scalar> im;
        };
        struct{
            data_arr data;
        };
    };

    quaternion(): quaternion(0, 0, 0, 0){};

    quaternion(scalar a, scalar b, scalar c, scalar d): a(a), b(b), c(c), d(d) {};

    quaternion(data_arr &v): a(v[0]), b(v[1]), c(v[2]), d(v[3]) {};

    quaternion(scalar re, vec3<scalar> im): re(re), im(im) {};

    quaternion(vec3<scalar> im): re(0), im(im) {};

    quaternion(std::vector<double> &v): a(v[0]), b(v[1]), c(v[2]), d(v[3]) {};

    //explicit quaternion(const data_arr &abcd) : data(abcd) {};

    quaternion &operator+=(const quaternion &rhs) {
        std::transform(data.begin(), data.end(), rhs.data.begin(), data.begin(), std::plus<scalar>());
        return *this;
    };

    quaternion &operator-=(const quaternion &rhs) {
        std::transform(data.begin(), data.end(), rhs.data.begin(), data.begin(), std::minus<scalar>());
        return *this;
    };

    quaternion operator*(const quaternion &other) const {
        return {
                re*other.re - im*other.im,
                re*other.im + other.re*im + im.cross(other.im)
        };
    };

    quaternion conj() {
        return {a, -b, -c, -d};
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    quaternion &operator+=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind1st(std::plus<scalar>(), a));
        return *this;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    quaternion &operator-=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind2nd(std::minus<scalar>(), a));
        return *this;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    quaternion &operator*=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind1st(std::multiplies<scalar>(), a));
        return *this;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    quaternion &operator/=(arithmetic a) {
        std::transform(data.begin(), data.end(), data.begin(), std::bind2nd(std::divides<scalar>(), a));
        return *this;
    };

    scalar norm() const {
        return std::sqrt(normSquared());
    };

    scalar normSquared() const {
        return std::inner_product(data.begin(), data.end(), data.begin(), static_cast<scalar>(0));
    };

    scalar operator[](std::size_t i) const {
        return data.at(i);
    };

    scalar &operator[](std::size_t i) {
        return data.at(i);
    };

    quaternion &invertElementWise() {
        for (auto i = 0; i < 3; ++i) {
            data[i] = static_cast<scalar>(1.) / data[i];
        }
        return *this;
    };

    bool operator==(const quaternion &rhs) const {
        return data == rhs.data;
    };

    bool operator!=(const quaternion &rhs) const {
        return data != rhs.data;
    };

    bool operator ==(const data_arr &rhs) const {
        return data == rhs;
    };

    bool operator ==(const std::vector<scalar> &rhs) const {
        return data == rhs;
    };
//    bool almostEquals(const quaternion &rhs) const {
//        bool result {true};
//        for(std::uint8_t i = 0; i < 3; ++i) {
//            const auto fp1 = fp::FloatingPoint<float>(data[i]);
//            const auto fp2 = fp::FloatingPoint<float>(rhs.data[i]);
//            auto x = fp::FloatingPoint<float>::DistanceBetweenSignAndMagnitudeNumbers(fp1.bits(), fp2.bits());
//            result &= fp1.AlmostEquals(fp2);
//        }
//        return result;
//    }

    friend std::ostream &operator<<(std::ostream &os, const quaternion &vec) {
        os << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ", " << vec[3] << ")";
        return os;
    };

    friend quaternion operator+(quaternion lhs, const quaternion &rhs) {
        lhs += rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend quaternion operator+(quaternion lhs, arithmetic rhs) {
        lhs += rhs;
        return lhs;
    };

    friend quaternion operator-(quaternion lhs, const quaternion &rhs) {
        lhs -= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend quaternion operator-(quaternion lhs, arithmetic rhs) {
        lhs -= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend quaternion operator/(quaternion lhs, arithmetic rhs) {
        lhs /= rhs;
        return lhs;
    };

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend quaternion operator*(quaternion lhs, const arithmetic rhs) {
        lhs *= rhs;
        return lhs;
    }

    template<typename arithmetic, typename detail::is_arithmetic_type<arithmetic> = 0>
    friend quaternion operator*(const arithmetic rhs, quaternion lhs) {
        lhs *= rhs;
        return lhs;
    }

};

// Other useful quaternion functions
template<typename scalar=double>
quaternion<scalar> angle2quaternion(const vec3<double> &phi) {
    double phinorm = phi.norm();
    if (phinorm != 0) {
        vec3<double> phiunit = phi / phinorm;
        double s = cos(0.5 * phinorm);
        double p = sin(0.5 * phinorm);
        return {s, p * phiunit[0], p * phiunit[1], p * phiunit[2]};
    } else {
        //returns unit quaternion (no rotation)
        return {1, 0, 0, 0};
    }
};