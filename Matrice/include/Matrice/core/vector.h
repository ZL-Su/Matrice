/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

#include <array>
#include "matrix.h"
#ifdef MATRICE_SIMD_ARCH
#include "arch/ixpacket.h"
#endif

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty, int _Dim = 2> 
class Vec_ : public Matrix_<_Ty, _Dim, compile_time_size<>::val_1>
{
	using _Myt = Vec_;
protected:
	using _Mybase = Matrix_<_Ty, _Dim, compile_time_size<>::val_1>;
	using const_initlist = typename _Mybase::const_initlist;
public:
	enum{rows_at_compiletime = _Dim, cols_at_compiletime = 1};
	using _Mybase::data;
	using typename _Mybase::value_t;
	using typename _Mybase::value_type;
	using reference = value_type&;
	using const_reference = const reference;

	MATRICE_GLOBAL_FINL Vec_() noexcept
		: _Mybase(0) {}
	MATRICE_GLOBAL_FINL Vec_(value_t _v) noexcept
		: _Mybase(_v) {}
	MATRICE_GLOBAL_FINL Vec_(value_t _x, value_t _y) noexcept
		: _Mybase({ _x, _y }) {}
	MATRICE_GLOBAL_FINL Vec_(const _Myt& _other) noexcept
		: _Mybase(_other) {}
	MATRICE_GLOBAL_FINL Vec_(const_initlist _list) noexcept
		: _Mybase(_list) {}
	MATRICE_GLOBAL_FINL Vec_(const _Mybase& _mtx) noexcept
		: _Mybase(_mtx) {}
	template<typename _Exp, MATRICE_ENABLE_IF(is_expression_v<_Exp>)>
	MATRICE_GLOBAL_FINL Vec_(const _Exp& _exp) noexcept {
		_exp.assign(*this); 
	}

	MATRICE_GLOBAL_FINL 
		reference operator[](size_t i) noexcept {
		return data()[i]; 
	}
	MATRICE_GLOBAL_FINL 
		const_reference operator[](size_t i)const noexcept {
		return data()[i]; 
	}
	MATRICE_GLOBAL_FINL 
		_Myt& operator= (const_initlist _list) noexcept {
		return static_cast<_Myt&>(_Mybase::operator= (_list)); 
	}
	MATRICE_GLOBAL_FINL 
		_Myt& operator= (const _Myt& _other) noexcept {
		return static_cast<_Myt&>(_Mybase::operator=(_other)); 
	}
	template<typename _Rval>
	MATRICE_GLOBAL_FINL 
		_Myt& operator= (const _Rval& _rval) noexcept {
		return static_cast<_Myt&>(_Mybase::operator= (_rval)); 
	}
	
	MATRICE_GLOBAL_FINL 
		operator typename _Mybase::pointer() noexcept {
		return data(); 
	}

	MATRICE_GLOBAL_FINL 
		_Myt& normalize(value_t _val = _Mybase::inf) noexcept {
		return static_cast<_Myt&>(
			_Mybase::operator = (_Mybase::normalize(_val))); 
	}

	MATRICE_GLOBAL_FINL
	/// <summary>
	/// \brief Compute dot-product with other vector. 
	/// </summary>
	/// <param name="'_other'">Other vector with the same type.</param>
	/// <returns></returns>
	value_t dot(const _Myt& _other) const noexcept {
		return _Mybase::dot(_other); 
	}

	MATRICE_GLOBAL_INL
	/// <summary>
	/// \brief Return evenly spaced values within a given interval. 
	/// The step is automatically according to the interval span and the vector dim.
	/// </summary>
	/// <param name= "'start'">Start of interval. The interval includes this value.</param>
	/// <param name= "'stop'">End of interval. The interval does not include this value.</param>
	/// <returns>Values are generated within the half-open interval [start, stop).</returns>
	static _Myt arange(value_t start, value_t end) noexcept {
		_Myt ret;
		const auto step = (end - start) / rows_at_compiletime;
		for (auto idx = 0; idx < rows_at_compiletime; ++idx) {
			ret[idx] = start + idx * step;
		}
		return ret;
	}
	
#ifdef _MSVC_LANG
	///<brief> properties </brief>
	__declspec(
		property(get = _x_getter, put = _x_setter)
		) reference x;
	MATRICE_GLOBAL_FINL 
		reference _x_getter()const noexcept { return data()[0]; }
	MATRICE_GLOBAL_FINL 
		void _x_setter(value_t _x)noexcept { data()[0] =_x; }
	__declspec(
		property(get = _y_getter, put = _y_setter)
		) reference y;
	MATRICE_GLOBAL_FINL 
		reference _y_getter()const noexcept { return data()[1]; }
	MATRICE_GLOBAL_FINL 
		void _y_setter(value_t _y)noexcept { data()[1] = _y; }
#else
	MATRICE_GLOBAL_INL const reference x() const noexcept {
		return data()[0];
	}
	MATRICE_GLOBAL_INL reference x() noexcept {
		return data()[0];
	}
	MATRICE_GLOBAL_INL const reference y() const noexcept {
		return data()[1];
	}
	MATRICE_GLOBAL_INL reference y() noexcept {
		return data()[1];
	}
#endif

	template<int rows>
	using lite = typename _Mybase::template 
		lite<rows, cols_at_compiletime>;
	template<int cols>
	using extend = typename _Mybase::template 
		lite<rows_at_compiletime, cols>;
};

template<typename _Ty> 
class Vec3_ MATRICE_NONHERITABLE : public Vec_<_Ty, 3>
{
	using _Myt = Vec3_;
	using _Mybase = Vec_<_Ty, 3>;
	using typename _Mybase::reference;
	using typename _Mybase::const_initlist;
public:
	using typename _Mybase::value_t;
	using typename _Mybase::value_type;
	using _Mybase::rows_at_compiletime;
	using _Mybase::cols_at_compiletime;
	using _Mybase::data;
	using _Mybase::x;
	using _Mybase::y;
	using _Mybase::Vec_;
	using _Mybase::operator=;
	using _Mybase::operator[];

	MATRICE_GLOBAL_FINL 
		Vec3_()noexcept {}
	MATRICE_GLOBAL_FINL
		Vec3_(value_t _x, value_t _y, value_t _z)noexcept
		: _Mybase({_x, _y, _z}) {}
	template<typename _Uy>
	MATRICE_GLOBAL_FINL 
		Vec3_(const Vec3_<_Uy>& _other)noexcept
		: Vec3_(_other.x, _other.y, _other.z) {}

	MATRICE_GLOBAL_FINL 
		_Myt& operator=(const _Myt& _other)noexcept { 
		return static_cast<_Myt&>(_Mybase::operator=(_other)); 
	}
	template<typename _Rval>
	MATRICE_GLOBAL_FINL 
		_Myt& operator=(const _Rval& _rval)noexcept { 
		return static_cast<_Myt&>(_Mybase::operator= (_rval)); 
	}

	MATRICE_GLOBAL_FINL 
		_Myt& normalize(value_t _val = _Mybase::inf)noexcept {
		return static_cast<_Myt&>(_Mybase::normalize(_val)); 
	}
	MATRICE_GLOBAL_FINL 
		value_t dot(const _Myt& _other)const noexcept { 
		return _Mybase::dot(_other); 
	}
	MATRICE_GLOBAL_FINL 
		_Myt cross(const _Myt& _rhs) const noexcept { 
		return _Myt(
			y*_rhs[2] - z*_rhs[1], 
			z*_rhs[0] - x * _rhs[2], 
			x*_rhs[1] - y * _rhs[0]
		);
	}

	/// <summary>
	/// \brief Span a 3-vector to a 3-by-3 skew-symmetric matrix.
	/// </summary>
	/// <returns> A 3-by-3 skew-symmetric matrix </returns>
	MATRICE_GLOBAL_INL 
		typename _Mybase::template extend<3>skew()const noexcept {
		return {0, -z, y, z, 0, -x, -y, x, 0};
	}

#ifdef _MSVC_LANG
	///<brief> properties </brief>
	__declspec(
		property(get = _z_getter, put = _z_setter)
		) reference z;
	MATRICE_GLOBAL_FINL 
		reference _z_getter()const noexcept { return data()[2]; }
	MATRICE_GLOBAL_FINL 
		void _z_setter(value_t _z)noexcept { data()[2] = _z; }
#else
	MATRICE_GLOBAL_INL const reference z() const noexcept {
		return data()[2];
}
	MATRICE_GLOBAL_INL reference z() noexcept {
		return data()[2];
	}
#endif

};

template<typename _Ty> 
class Vec4_ MATRICE_NONHERITABLE : public Vec_<_Ty, 4>
{
	using _Myt = Vec4_;
	using _Mybase = Vec_<_Ty, 4>;
	using typename _Mybase::reference;
	using typename _Mybase::const_initlist;
public:
	using typename _Mybase::value_t;
	using typename _Mybase::value_type;
	using _Mybase::rows_at_compiletime;
	using _Mybase::cols_at_compiletime;
	using _Mybase::data;
	using _Mybase::x;
	using _Mybase::y;
	using _Mybase::Vec_;
	using _Mybase::operator=;
	using _Mybase::operator[];

	MATRICE_GLOBAL_FINL 
		Vec4_(value_t _x, value_t _y, value_t _z) noexcept
		: _Mybase({ _x, _y, _z, 1}) {}
	MATRICE_GLOBAL_FINL 
		Vec4_(value_t _x,value_t _y,value_t _z,value_t _w) noexcept
		: _Mybase({ _x, _y, _z, _w }) {}
	template<typename _Uy>
	MATRICE_GLOBAL_FINL 
		Vec4_(const Vec4_<_Uy>& _other) noexcept
		: Vec4_(_other.x, _other.y, _other.z, _other.w) {}

	template<typename _Rval>
	MATRICE_GLOBAL_FINL 
		_Myt& operator= (const _Rval& _rval) noexcept { 
		return static_cast<_Myt&>(_Mybase::operator= (_rval)); 
	}

	MATRICE_GLOBAL_FINL 
		_Myt& normalize(value_t _val = _Mybase::inf) noexcept {
		return static_cast<_Myt&>(_Mybase::normalize(_val)); 
	}

#ifdef _MSVC_LANG
	///<brief> properties </brief>
	__declspec(
		property(get = _z_getter, put = _z_setter)
		) reference z;
	MATRICE_GLOBAL_FINL 
		reference _z_getter()const noexcept { return data()[2]; }
	MATRICE_GLOBAL_FINL 
		void _z_setter(value_t _z) noexcept { data()[2] = _z; }
	__declspec(
		property(get = _w_getter, put = _w_setter)
		) reference w;
	MATRICE_GLOBAL_FINL 
		reference _w_getter()const noexcept { return data()[3]; }
	MATRICE_GLOBAL_FINL 
		void _w_setter(value_t _w) noexcept { data()[3] = _w; }
#else
	MATRICE_GLOBAL_INL const reference z() const noexcept {
		return data()[2];
	}
	MATRICE_GLOBAL_INL reference z() noexcept {
		return data()[2];
	}
	MATRICE_GLOBAL_INL const reference w() const noexcept {
		return data()[3];
	}
	MATRICE_GLOBAL_INL reference w() noexcept {
		return data()[3];
	}
#endif
};
_DETAIL_END

// \brief Generic managed vector type
template<typename _Ty, int _Dim, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using Vec_ = detail::Vec_<_Ty, _Dim>;
template<typename _Ty, int _Dim>
struct is_fxdvector<Vec_<_Ty, _Dim>> : std::true_type {};

// \brief Managed vector type with 2 entities: x, y
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using Vec2_ = detail::Vec_<_Ty, 2>;
template<typename _Ty>
struct is_fxdvector<Vec2_<_Ty>> : std::true_type {};

// \brief Managed vector type with 3 entities: x, y, z
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using Vec3_ = detail::Vec3_<_Ty>;
template<typename _Ty>
struct is_fxdvector<Vec3_<_Ty>> : std::true_type {};

// \brief Managed vector type with 4 entities: x, y, z, w
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using Vec4_ = detail::Vec4_<_Ty>;
template<typename _Ty>
struct is_fxdvector<Vec4_<_Ty>> : std::true_type {};

// \brief Dispatch vector type according to a given dimensionality _Dim.
template<typename _Ty, size_t _Dim> struct auto_vector {
	using type = Vec_<_Ty, _Dim>;
};
template<typename _Ty> struct auto_vector<_Ty, 2> { 
	using type = Vec2_<_Ty>; 
};
template<typename _Ty> struct auto_vector<_Ty, 3> { 
	using type = Vec3_<_Ty>; 
};
template<typename _Ty> struct auto_vector<_Ty, 4> { 
	using type = Vec4_<_Ty>; 
};
template<typename _Ty> struct auto_vector<_Ty, ::dynamic> { 
	using type = std::vector<_Ty>; 
};

/// <summary>
/// ALIAS, auto_vector type
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty, size_t _Dim>
using auto_vector_t = typename auto_vector<_Ty, _Dim>::type;

/**
 *\brief Matrice also regards std::array as a fixed vector type.
 */
template<typename _Ty, int _Dim>
struct is_fxdvector<std::array<_Ty,_Dim>> : std::true_type {};

template<size_t _Dim, typename _Ty> MATRICE_GLOBAL_INL 
auto cross_prod_matrix(const Vec3_<_Ty>& v) noexcept {
	static_assert(_Dim != 3 || _Dim != 4,
		"The _Dim in function cross_prod_matrix<_Dim, _Ty> must be 3 or 4.");
	using _Rety = detail::Matrix_<_Ty, _Dim, _Dim>;

	MATRICE_CONSTEXPR_IF(_Dim == 3)
		return _Rety{ 0, -v.z, v.y, v.z, 0, -v.x, -v.y, v.x, 0 };
	else
		return _Rety{ 0, -v.z, v.y, 0, v.z, 0, -v.x, 0, -v.y, v.x, 0, 0, 0, 0, 0, 0 };
}

/**
 * \brief Concatenate two vectors as a matrix in column-by-column order.
 */
template<typename _Ty, size_t _N> MATRICE_GLOBAL_INL
auto concat(const Vec_<_Ty, _N>& prev, const Vec_<_Ty, _N>& next) noexcept {
	typename Vec_<_Ty, _N>::template extend<2> _Ret;
	_Ret.cview(0) = prev, _Ret.cview(1) = next;
	return _Ret;
}

/**
 * \brief Span a 3-vector to a skew-symmetric matrix.
 * \param "v" the input 3d vector.
 * \returns A matrix with type of Matrix_<_Ty, 3, 3>.
 */
template<typename _Ty> MATRICE_GLOBAL_INL
auto skew(const Vec3_<_Ty>& v) noexcept {
	return v.skew();
}
template<typename _Ty> MATRICE_GLOBAL_INL
auto skew(const initlist<_Ty> v) noexcept {
	return Vec3_<_Ty>(v).skew();
}

DGE_MATRICE_END