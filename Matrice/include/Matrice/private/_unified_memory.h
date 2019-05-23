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
#include <type_traits>
#include "../util/_macros.h"
#include "_type_traits.h"
#include "_memory.h"

MATRICE_PRIVATE_BEGIN
template<typename T> struct memory_base_v1 {
	using value_t = T;
	using pointer = value_t*;
	using const_pointer = const pointer;
};
///<!-- generic class for device memory operation-->
template<typename _Scalar, Location _From, Location _To, size_t _Opt>
struct unified_sync : memory_base_v1<_Scalar> {
	enum { From = _From, To = _To, Option = _Opt };
	static_assert(Option == LINEAR + COPY, "Only for host linear memory copy!");
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	MATRICE_GLOBAL_FINL static pointer op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols = 1, size_t _1 = 1)
	{ return fill_mem(_Src, _Dst, _Rows*_Cols); }
};
template<typename _Scalar, Location _Host1, Location _Host2>
struct unified_sync<_Scalar, _Host1, _Host2, LINEAR + MOVE>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = _Host1, To = _Host2, Option = LINEAR + MOVE };
	MATRICE_HOST_FINL static pointer op(pointer _Dst, pointer _Src, size_t _1 = 1, size_t _2 = 1, size_t _3 = 1)
	{_Dst = _Src; _Src = nullptr; return (_Dst); }
};
template<typename _Scalar, Location _Host1, Location _Host2>
struct unified_sync<_Scalar, _Host1, _Host2, LINEAR + SHARED>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = _Host1, To = _Host2, Option = LINEAR + SHARED};
	MATRICE_HOST_FINL static pointer op(pointer _Dst, pointer _Src, size_t _1 = 1, size_t _2 = 1, size_t _3 = 1)
	{ return (_Dst = _Src); }
};
template<typename _Scalar>
struct unified_sync<_Scalar, OnStack, OnStack, LINEAR + SHARED>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = OnStack, To = OnStack, Option = LINEAR + SHARED };
	MATRICE_GLOBAL_FINL static pointer op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols = 1, size_t _1 = 1)
	{ return fill_mem(_Src, _Dst, _Rows*_Cols); }
};
template<typename _Scalar, Location _Host>
struct unified_sync<_Scalar, _Host, OnDevice, LINEAR>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = _Host, To = OnDevice, Option = LINEAR };
	MATRICE_GLOBAL static pointer op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols = 1, size_t _1 = 1);
};
template<typename _Scalar, Location _Host>
struct unified_sync<_Scalar, _Host, OnDevice, PITCHED>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = _Host, To = OnDevice, Option = PITCHED };
	MATRICE_GLOBAL static pointer op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _Pytes);
};
template<typename _Scalar, Location _Host>
struct unified_sync<_Scalar, OnDevice, _Host, LINEAR>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = OnDevice, To = _Host, Option = LINEAR };
	MATRICE_GLOBAL static pointer op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols = 1, size_t _1 = 1);
};
template<typename _Scalar, Location _Host>
struct unified_sync<_Scalar, OnDevice, _Host, PITCHED>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = OnDevice, To = _Host, Option = PITCHED };
	MATRICE_GLOBAL static pointer op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _Pytes);
};
template<typename _Scalar>
struct unified_sync<_Scalar, OnDevice, OnDevice, LINEAR>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = OnDevice, To = OnDevice, Option = LINEAR };
	MATRICE_GLOBAL static pointer op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols = 1, size_t _1 = 1);
};
template<typename _Scalar>
struct unified_sync<_Scalar, OnDevice, OnDevice, PITCHED>
	: memory_base_v1<_Scalar>
{
	using typename memory_base_v1<_Scalar>::pointer;
	using typename memory_base_v1<_Scalar>::const_pointer;
	enum { From = OnDevice, To = OnDevice, Option = PITCHED };
	MATRICE_GLOBAL static pointer op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _Pytes);
};

///<brief> unified memory release </brief>
template<typename _Ty, Location _Loc> struct unified_free {};
template<typename _Ty> struct unified_free<_Ty, OnHeap> {
	MATRICE_GLOBAL_FINL void operator() (_Ty* _Data) { aligned_free(_Data); }
};
template<typename _Ty> struct unified_free<_Ty, OnDevice> {
	MATRICE_GLOBAL_FINL void operator() (_Ty* _Data) { device_free(_Data); }
};
template<typename _Ty> struct unified_free<_Ty, OnGlobal> {
	MATRICE_GLOBAL_FINL void operator() (_Ty* _Data) { device_free(_Data); }
};

///<brief> unified memory set </brief>
template<typename _Ty, size_t _Opt> struct unified_fill {
	using FwdIt = _Ty * ;
	MATRICE_HOST_FINL static void op(FwdIt _First, const _Ty& _Val, size_t _Size, size_t _1 = 1, size_t _2 = 1) {
		std::fill(_First, _First + _Size*_1, _Val);
	}
};
template<typename _Ty> struct unified_fill<_Ty, OnDevice + LINEAR> {
	using FwdIt = _Ty * ;
	MATRICE_GLOBAL_FINL static void op(FwdIt _First, const _Ty& _Val, size_t _Size, size_t _1 = 1, size_t _2 = 1) {
		//cudaMemset(_First, _Val, _Size*_1*type_bytes<_Ty>::value);
		std::runtime_error("Oops, filling linear memory on device has not been implemented.");
	}
};
template<typename _Ty> struct unified_fill<_Ty, OnDevice + PITCHED> {
	using FwdIt = _Ty * ;
	MATRICE_GLOBAL_FINL static void op(FwdIt _First, const _Ty& _Val, size_t _Rows, size_t _Cols, size_t _Pitch) {
		//cudaMemset2D(_First, _Pitch, _Val, _Cols*type_bytes<_Ty>::value, _Rows);
		std::runtime_error("Oops, filling pitched device memory has not been implemented.");
	}
};
MATRICE_PRIVATE_END