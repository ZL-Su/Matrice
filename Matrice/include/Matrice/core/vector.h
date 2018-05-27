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

#include "matrix.h"
#ifndef __vec_global_inl
#define __vec_global_inl MATRICE_GLOBAL_INL
#endif

MATRICE_NAMESPACE_BEGIN_TYPES

template<typename _Ty, int _Dim = 2> 
class Vec_ : public Matrix_<_Ty, _Dim, compile_time_size<>::val_1>
{
	using const_my = const Vec_;
	using const_my_ref = const_my&;
protected:
	using _Base = Matrix_<_Ty, _Dim, compile_time_size<>::val_1>;
	using const_init_list = typename _Base::const_init_list;
	typename _Base::pointer _Data = _Base::data();
public:
	using value_t = typename _Base::value_t;
	using const_value = const value_t;
	using reference = value_t & ;
	using const_reference = const reference;
	__vec_global_inl Vec_() : _Base({ 0 }) {}
	__vec_global_inl Vec_(const_value _v) : _Base({ _v }) {}
	__vec_global_inl Vec_(const_value _x, const_value _y) : _Base({ _x, _y }) {}
	__vec_global_inl Vec_(const_my_ref _other) : _Base(_other) {}
	__vec_global_inl Vec_(const_init_list _list) : _Base(_list) {}
	template<typename _Expr>
	__vec_global_inl Vec_(const _Expr& _xpr) : _Base(_xpr) {}

	__vec_global_inl Vec_& operator= (const_init_list _list)
	{ return static_cast<Vec_&>(_Base::operator= (_list)); }
	__vec_global_inl Vec_& operator= (const_my_ref _other)
	{ return static_cast<Vec_&>(_Base::operator=(_other)); }
	template<typename _Rval>
	__vec_global_inl Vec_& operator= (const _Rval& _rval)
	{ return static_cast<Vec_&>(_Base::operator= (_rval)); }
	__vec_global_inl Vec_& normalize(const value_t _val = _Base::inf)
	{ return static_cast<Vec_&>(_Base::operator = (_Base::normalize(_val))); }
	__vec_global_inl value_t dot(const_my_ref _other) const
	{ return _Base::dot(_other); }
	
	///<brief> properties </brief>
	__declspec(property(get = _x_getter, put = _x_setter)) reference x;
	__vec_global_inl reference _x_getter() const { return _Data[0]; }
	__vec_global_inl void _x_setter(value_t _x) { _Data[0] =_x; }
	__declspec(property(get = _y_getter, put = _y_setter)) reference y;
	__vec_global_inl reference _y_getter() const { return _Data[1]; }
	__vec_global_inl void _y_setter(value_t _y) { _Data[1] = _y; }

};
template<typename _Ty> class Vec3_ final : public Vec_<_Ty, 3>
{
	using _Base = Vec_<_Ty, 3>;
	using const_my_ref = const Vec3_&;
	using typename _Base::const_value;
	using typename _Base::reference;
	using typename _Base::const_init_list;
	using _Base::_Data;
public:
	using typename _Base::value_t;
	using _Base::operator=;
	using _Base::x;
	using _Base::y;
	__vec_global_inl Vec3_() : _Base() {}
	__vec_global_inl Vec3_(const_value _v) : _Base(_v) {}
	__vec_global_inl Vec3_(const_value _x, const_value _y) : _Base({ _x, _y, value_t(1) }) {}
	__vec_global_inl Vec3_(const_value _x, const_value _y, const_value _z) : _Base({_x, _y, _z}) {}
	__vec_global_inl Vec3_(const_init_list _list) : _Base(_list) {}
	template<typename _Type>
	__vec_global_inl Vec3_(const _Type& _other) : _Base(_other) {}

	__vec_global_inl Vec3_& operator= (const_my_ref _other)
	{ return static_cast<Vec3_&>(_Base::operator=(_other)); }
	template<typename _Rval>
	__vec_global_inl Vec3_& operator= (const _Rval& _rval)
	{ return static_cast<Vec3_&>(_Base::operator= (_rval)); }
	__vec_global_inl Vec3_& normalize(const value_t _val = _Base::inf)
	{ return static_cast<Vec3_&>(_Base::normalize(_val)); }
	__vec_global_inl value_t dot(const_my_ref _other) const
	{ return _Base::dot(_other); }
	__vec_global_inl Vec3_ cross(const_my_ref _rhs) const
	{ return Vec3_(y*_rhs.z - z*_rhs.y, z*_rhs.x - x * _rhs.z, x*_rhs.y - y * _rhs.x); }
	__declspec(property(get = _z_getter, put = _z_setter)) reference z;
	__vec_global_inl reference _z_getter() const { return _Data[2]; }
	__vec_global_inl void _z_setter(value_t _z) { _Data[2] = _z; }
};
template<typename _Ty> class Vec4_ final : public Vec_<_Ty, 4>
{
	using _Base = Vec_<_Ty, 4>;
	using typename _Base::const_value;
	using typename _Base::reference;
	using typename _Base::const_init_list;
	using _Base::_Data;
public:
	using typename _Base::value_t;
	using _Base::operator=;
	using _Base::x;
	using _Base::y;
	__vec_global_inl Vec4_() : _Base() {}
	__vec_global_inl Vec4_(const_value _v) : _Base(_v) {}
	__vec_global_inl Vec4_(const_value _x, const_value _y, const_value _z) : _Base({ _x, _y, _z, 1}) {}
	__vec_global_inl Vec4_(const_value _x, const_value _y, const_value _z, const_value _w) : _Base({ _x, _y, _z, _w }) {}
	__vec_global_inl Vec4_(const_init_list _list) : _Base(_list) {}
	template<typename _Expr>
	__vec_global_inl Vec4_(const _Expr& _xpr) : _Base(_xpr) {}
	template<typename _Rval>
	__vec_global_inl Vec4_& operator= (const _Rval& _rval)
	{ return static_cast<Vec4_&>(_Base::operator= (_rval)); }
	__vec_global_inl Vec4_& normalize(const value_t _val = _Base::inf)
	{ return static_cast<Vec4_&>(_Base::normalize(_val)); }
	__declspec(property(get = _z_getter, put = _z_setter)) reference z;
	__vec_global_inl reference _z_getter() const { return _Data[2]; }
	__vec_global_inl void _z_setter(value_t _z) { _Data[2] = _z; }
	__declspec(property(get = _w_getter, put = _w_setter)) reference w;
	__vec_global_inl reference _w_getter() const { return _Data[3]; }
	__vec_global_inl void _w_setter(value_t _w) { _Data[3] = _w; }
};
MATRICE_NAMESPACE_END_TYPES