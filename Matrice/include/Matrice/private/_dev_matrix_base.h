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
#include "../util/utils.h"
#if (defined MATRICE_ENABLE_CUDA && !defined __disable_cuda__)
#include "_decl_dev_funcs.h"
#include "_devops.h"
#include "_unified_memory.h"

MATRICE_NAMESPACE_BEGIN_TYPES
template<typename _Ty, int _M, int _N> class Matrix_;
MATRICE_NAMESPACE_END_TYPES

MATRICE_DEVICE_BEGIN

// \CLASS TEMPLATE to make base class for device matrix
template<typename _Ty> class Base_
{
	using _Myt = Base_;
	using value_t = _Ty;
	using _Mydt = types::Matrix_<value_t, -1, -1>;
	using int_ptr_t = std::add_pointer_t<int>;
	using pointer = std::add_pointer_t<value_t>;
public:
	static constexpr size_t threads_per_block = 128;

	MATRICE_GLOBAL_INL Base_() = default;
	MATRICE_GLOBAL_INL Base_(pointer pdev, size_t* p, int_ptr_t w, int_ptr_t h)
		:_Ptr(pdev), _P(p), _W(w), _H(h) {};

	//<brief> synchronize all device threads </brief>
	MATRICE_HOST_INL void sync() { 
		_Future = std::async(std::launch::async,[&]{_Sync_impl();}); 
	}

	//<brief> retrieve data from device memory </brief>
	MATRICE_HOST_INL void fetch(pointer dst) { 
		if (_Future.valid()) _Future.get(); _Dnload_impl(dst); 
	}
	template<class _Arg, MATRICE_ENABLE_IF(is_class_v<_Arg>)>
	MATRICE_HOST_INL void fetch(_Arg& dst) { 
		if (_Future.valid()) _Future.get(); _Dnload_impl(dst.data()); 
	}
	/**
	 * \sync data from device to host deglom::Matrix_ type
 	 */
	template<int _M = 0, int _N = _M>
	MATRICE_HOST_INL auto fetch() {
		using _Rety = types::Matrix_<value_t, _M, _N>;
		_Rety _Ret(*_H, *_W);

		if (_Future.valid()) _Future.get();
		_Dnload_impl(_Ret.data());

		return std::forward<_Rety>(_Ret);
	}

	//<brief> copy data from host to device memory </brief>
	MATRICE_GLOBAL_INL pointer operator= (const pointer _src) {
		_Upload_impl(_src); return (_Ptr);
	}
	template<typename _Rhs, MATRICE_ENABLE_IF(is_class_v<_Rhs>)>
	MATRICE_GLOBAL_INL _Mydt& operator= (const _Rhs& _src) {
		_Upload_impl(_src.data()); return *static_cast<_Mydt*>(this);
	}

	/**
	 * if the device has the compute_35 architecture or above, using MATRICE_GLOBAL_INL
	 */
	MATRICE_HOST_INL pointer operator+(const _Myt& _other);
	MATRICE_HOST_INL pointer operator-(const _Myt& _other);
	MATRICE_HOST_INL pointer operator*(const _Myt& _other);
	MATRICE_HOST_INL pointer operator/(const _Myt& _other);

	MATRICE_HOST_INL value_t reduce();
	
	friend _Mydt operator+(const value_t _Left, const _Mydt& _Right);
	friend _Mydt operator-(const value_t _Left, const _Mydt& _Right);
	friend _Mydt operator*(const value_t _Left, const _Mydt& _Right);
	friend _Mydt operator/(const value_t _Left, const _Mydt& _Right);

protected:
	MATRICE_HOST_INL void _Sync_impl() { privt::_Device_sync<0>(); }

private:
	template<typename... Args>
	MATRICE_HOST_INL void _Upload_impl(pointer args) { device_memcpy<value_t, 1>::impl(args, _Ptr, size_t(*_W), size_t(*_H), *_P); }
	template<typename... Args>
	MATRICE_HOST_INL void _Dnload_impl(pointer args) { device_memcpy<value_t, 2>::impl(args, _Ptr, size_t(*_W), size_t(*_H), *_P); }

	std::future<void> _Future; 
	pointer _Ptr;
	size_t* _P = nullptr; 
	int_ptr_t  _W = nullptr, _H = nullptr;
};
MATRICE_DEVICE_END

#endif