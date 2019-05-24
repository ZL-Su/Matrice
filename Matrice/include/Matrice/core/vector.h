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

#ifndef vec_global_inl
#define vec_global_inl constexpr MATRICE_GLOBAL_FINL
#endif

MATRICE_NAMESPACE_BEGIN_TYPES
template<typename _Ty, int _Dim = 2> 
class Vec_ : public Matrix_<_Ty, _Dim, compile_time_size<>::val_1>
{
	using const_my = const Vec_;
	using my_const_ref = const_my&;
protected:
	using _Base = Matrix_<_Ty, _Dim, compile_time_size<>::val_1>;
	using const_initlist = typename _Base::const_initlist;
public:
	enum{CompileTimeRows = _Dim, CompileTimeCols = 1};
	using _Base::data;
	using typename _Base::value_t;
	using typename _Base::value_type;
	using const_value = const value_t;
	using reference = value_t & ;
	using const_reference = const reference;

	vec_global_inl Vec_() : _Base({ 0 }) {}
	vec_global_inl Vec_(const_value _v) : _Base({ _v }) {}
	vec_global_inl Vec_(const_value _x, const_value _y) : _Base({ _x, _y }) {}
	vec_global_inl Vec_(my_const_ref _other) : _Base(_other) {}
	vec_global_inl Vec_(const_initlist _list) : _Base(_list) {}
	template<typename _Exp, MATRICE_ENABLE_IF(is_expression_v<_Exp>)>
	vec_global_inl Vec_(const _Exp& _exp) { _exp.assign(*this); }

	vec_global_inl reference operator[] (size_t i) { return data()[i]; }
	vec_global_inl const_reference operator[](size_t i)const { return data()[i]; }
	vec_global_inl Vec_& operator= (const_initlist _list)
	{ return static_cast<Vec_&>(_Base::operator= (_list)); }
	vec_global_inl Vec_& operator= (my_const_ref _other)
	{ return static_cast<Vec_&>(_Base::operator=(_other)); }
	template<typename _Rval>
	vec_global_inl Vec_& operator= (const _Rval& _rval)
	{ return static_cast<Vec_&>(_Base::operator= (_rval)); }
	
	vec_global_inl operator typename _Base::pointer() { return data(); }

	vec_global_inl Vec_& normalize(const value_t _val = _Base::inf)
	{ return static_cast<Vec_&>(_Base::operator = (_Base::normalize(_val))); }
	vec_global_inl value_t dot(my_const_ref _other) const
	{ return _Base::dot(_other); }
	
	///<brief> properties </brief>
	__declspec(property(get = _x_getter, put = _x_setter)) reference x;
	vec_global_inl reference _x_getter() const { return data()[0]; }
	vec_global_inl void _x_setter(value_t _x) { data()[0] =_x; }
	__declspec(property(get = _y_getter, put = _y_setter)) reference y;
	vec_global_inl reference _y_getter() const { return data()[1]; }
	vec_global_inl void _y_setter(value_t _y) { data()[1] = _y; }

};
template<typename _Ty> class Vec3_ final : public Vec_<_Ty, 3>
{
	using _Base = Vec_<_Ty, 3>;
	using my_const_ref = const Vec3_&;
	using typename _Base::const_value;
	using typename _Base::reference;
	using typename _Base::const_initlist;
public:
	using _Base::CompileTimeRows;
	using _Base::CompileTimeCols;
	using typename _Base::value_t;
	using typename _Base::value_type;
	using _Base::operator=;
	using _Base::operator[];
	using _Base::data;
	using _Base::x;
	using _Base::y;
	using _Base::Vec_;

	vec_global_inl Vec3_(const_value _x, const_value _y, const_value _z) 
		: _Base({_x, _y, _z}) {}
	template<typename _Uy>
	vec_global_inl Vec3_(const Vec3_<_Uy>& _other) 
		: Vec3_(_other.x, _other.y, _other.z) {}

	vec_global_inl Vec3_& operator= (my_const_ref _other)
	{ return static_cast<Vec3_&>(_Base::operator=(_other)); }
	template<typename _Rval>
	vec_global_inl Vec3_& operator= (const _Rval& _rval)
	{ return static_cast<Vec3_&>(_Base::operator= (_rval)); }

	vec_global_inl Vec3_& normalize(const value_t _val = _Base::inf)
	{ return static_cast<Vec3_&>(_Base::normalize(_val)); }
	vec_global_inl value_t dot(my_const_ref _other) const
	{ return _Base::dot(_other); }
	vec_global_inl Vec3_ cross(my_const_ref _rhs) const
	{ return Vec3_(y*_rhs[2] - z*_rhs[1], z*_rhs[0] - x * _rhs[2], x*_rhs[1] - y * _rhs[0]); }
	__declspec(property(get = _z_getter, put = _z_setter)) reference z;
	vec_global_inl reference _z_getter() const { return data()[2]; }
	vec_global_inl void _z_setter(value_t _z) { data()[2] = _z; }
};

template<typename _Ty> class Vec4_ final : public Vec_<_Ty, 4>
{
	using _Base = Vec_<_Ty, 4>;
	using typename _Base::const_value;
	using typename _Base::reference;
	using typename _Base::const_initlist;
public:
	using _Base::CompileTimeRows;
	using _Base::CompileTimeCols;
	using typename _Base::value_t;
	using typename _Base::value_type;
	using _Base::operator=;
	using _Base::operator[];
	using _Base::data;
	using _Base::x;
	using _Base::y;
	using _Base::Vec_;

	vec_global_inl Vec4_(const_value _x, const_value _y, const_value _z) 
		: _Base({ _x, _y, _z, 1}) {}
	vec_global_inl Vec4_(const_value _x,const_value _y,const_value _z,const_value _w) 
		: _Base({ _x, _y, _z, _w }) {}
	template<typename _Uy>
	vec_global_inl Vec4_(const Vec4_<_Uy>& _other)
		: Vec4_(_other.x, _other.y, _other.z, _other.w) {}

	template<typename _Rval>
	vec_global_inl Vec4_& operator= (const _Rval& _rval)
	{ return static_cast<Vec4_&>(_Base::operator= (_rval)); }

	vec_global_inl Vec4_& normalize(const value_t _val = _Base::inf)
	{ return static_cast<Vec4_&>(_Base::normalize(_val)); }

	__declspec(property(get = _z_getter, put = _z_setter)) reference z;
	vec_global_inl reference _z_getter() const { return data()[2]; }
	vec_global_inl void _z_setter(value_t _z) { data()[2] = _z; }
	__declspec(property(get = _w_getter, put = _w_setter)) reference w;
	vec_global_inl reference _w_getter() const { return data()[3]; }
	vec_global_inl void _w_setter(value_t _w) { data()[3] = _w; }
};

MATRICE_NAMESPACE_END_TYPES

DGE_MATRICE_BEGIN
// \a generic managed vector type
template<typename _Ty, int _Dim,
	typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
using Vec_ = types::Vec_<_Ty, _Dim>;
template<typename _Ty, int _Dim>
struct is_fxdvector<Vec_<_Ty, _Dim>> : std::true_type {};

// \managed vector type with 2 entities: x, y
template<typename _Ty, 
	typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
using Vec2_ = types::Vec_<_Ty, 2>;
template<typename _Ty>
struct is_fxdvector<Vec2_<_Ty>> : std::true_type {};

// \managed vector type with 3 entities: x, y, z
template<typename _Ty,
	typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
using Vec3_ = types::Vec3_<_Ty>;
template<typename _Ty>
struct is_fxdvector<Vec3_<_Ty>> : std::true_type {};

// \managed vector type with 4 entities: x, y, z, w
template<typename _Ty,
	typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
using Vec4_ = types::Vec4_<_Ty>;
template<typename _Ty>
struct is_fxdvector<Vec4_<_Ty>> : std::true_type {};

/**
 *\brief Matrice regards the std::array is also a fixed vector type
 */
template<typename _Ty, int _Dim>
struct is_fxdvector<std::array<_Ty,_Dim>> : std::true_type {};

template<size_t _Dim, typename _Ty>
vec_global_inl types::Matrix_<_Ty, _Dim, _Dim> cross_prod_matrix(const Vec3_<_Ty>& v) {
	using return_t = types::Matrix_<_Ty, _Dim, _Dim>;
	MATRICE_CONSTEXPR_IF(_Dim == 3)
	return return_t{ 0, -v.z, v.y, v.z, 0, -v.x, -v.y, v.x, 0 };
	else
	return return_t{ 0, -v.z, v.y, 0, v.z, 0, -v.x, 0, -v.y, v.x, 0, 0, 0, 0, 0, 0 };
	static_assert(_Dim != 3 || _Dim != 4, "The _Dim in function cross_prod_matrix<_Dim> must be 3 or 4.");
}
DGE_MATRICE_END