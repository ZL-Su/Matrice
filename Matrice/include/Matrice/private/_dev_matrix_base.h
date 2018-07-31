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
#include <future>
#include "../util/_macros.h"
#if (defined __enable_cuda__ && !defined __disable_cuda__)
#include "_decl_dev_funcs.h"
#include "_devops.h"
#include "_unified_memory.h"

MATRICE_NAMESPACE_BEGIN_TYPES
template<typename _Ty, int _M, int _N> class Matrix_;
MATRICE_NAMESPACE_END_TYPES

MATRICE_DEVICE_BEGIN
template<typename _Ty, typename _Derived = types::Matrix_<_Ty, -1, -1>> 
class Base_
{
	using size_t = std::size_t;
	using intp_t = std::add_pointer_t<int>;
	using value_t = _Ty;
	using pointer = std::add_pointer_t<value_t>;
	using derived_t = _Derived;
public:
	MATRICE_GLOBAL_INL Base_() = default;
	MATRICE_GLOBAL_INL Base_(pointer pdev, size_t* p, intp_t w, intp_t h):_Ptr(pdev), _P(p), _W(w), _H(h) {};

	//<brief> synchronize all device threads </brief>
	MATRICE_HOST_INL void sync() { _Future = std::async(std::launch::async, [&] { _Sync_impl(); }); }
	//<brief> retrieve data from device memory </brief>
	MATRICE_HOST_INL void fetch(pointer dst) { if (_Future.valid()) _Future.get(); _Dnload_impl(dst); }
	template<typename _Arg, typename = std::enable_if_t<std::is_class_v<_Arg>>>
	MATRICE_HOST_INL void fetch(_Arg& dst){ if (_Future.valid()) _Future.get(); _Dnload_impl(dst.data()); }
	//<brief> copy data from host to device memory </brief>
	MATRICE_GLOBAL_INL _Derived& operator= (const pointer _src) { _Upload_impl(_src); return *static_cast<_Derived*>(this); }
	template<typename _Arg, typename = std::enable_if_t<std::is_class_v<_Arg>>>
	MATRICE_GLOBAL_INL _Derived& operator= (const _Arg& _src) { _Upload_impl(_src.data()); return *static_cast<_Derived*>(this); }

	MATRICE_GLOBAL_INL _Derived& operator+(const _Derived& _other);
private:
	MATRICE_HOST_INL void _Sync_impl() { privt::_Device_sync<0>(); }
	template<typename... Args>
	MATRICE_HOST_INL void _Upload_impl(pointer args) { device_memcpy<_Ty, 1>::impl(args, _Ptr, size_t(*_W), size_t(*_H), *_P); }
	template<typename... Args>
	MATRICE_HOST_INL void _Dnload_impl(pointer args) { device_memcpy<_Ty, 2>::impl(args, _Ptr, size_t(*_W), size_t(*_H), *_P); }

	std::future<void> _Future; 
	pointer _Ptr;
	size_t* _P = nullptr; 
	intp_t  _W = nullptr, _H = nullptr;
};
MATRICE_DEVICE_END

#endif