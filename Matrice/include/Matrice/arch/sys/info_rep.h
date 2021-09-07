/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/
#pragma once
#include "util/_macros.h"
#include <stdio.h>

DGE_MATRICE_BEGIN
namespace sys {
using byte_pointer = unsigned char*;
MATRICE_HOST_INL void _Show_bytes(byte_pointer _Begin, size_t _Len) noexcept {
	for (auto _Idx = 0; _Idx < _Len; ++_Idx) {
		printf(" %.2x", _Begin[_Idx]);
	}
	printf("\n");
}

/// <summary>
/// \brief Print the byte representation of program object '_Val'.
/// </summary>
template<typename _Ty>
MATRICE_HOST_INL void byte_rep(_Ty _Val) noexcept {
	_Show_bytes(byte_pointer(&_Val), sizeof(_Ty));
}

/// <summary>
/// \brief Swap to the values stored at the locations denoted by pointers x and y.
/// </summary>
template<typename _Ty>
MATRICE_HOST_INL void inplace_swap(_Ty* _Left, _Ty* _Right) noexcept {
	*_Right = *_Left ^ *_Right;
	*_Left = *_Left ^ *_Right;
	*_Right = *_Left ^ *_Right;
}

/// <summary>
/// \brief Reverse element order in a given array '_Array'.
/// </summary>
template<typename _Ty, size_t _N>
MATRICE_HOST_INL void reverse_array(_Ty(&_Array)[_N]) noexcept {
	auto _First = 0ull, _Last = _N - 1;
	for (; _First < _Last; ++_First, --_Last) {
		inplace_swap(&_Array[_First], &_Array[_Last]);
	}
}

/// <summary>
/// \brief Convert arrays of zeros and ones of length _W to nonnegative integers.
/// </summary>
/// <param name="_Bits">: Array of zeros and ones.</param>
/// <returns> A nonnegative integer. </returns>
template<size_t _W>
MATRICE_HOST_INL uint32_t B2U(bool(&_Bits)[_W]) noexcept {
	auto _Ret = uint32_t(0);
	for (auto _Idx = 0; _Idx < _W; ++_Idx) {
		_Ret += _Bits[_Idx] * pow(2, _W - 1 - _Idx);
	}
	return _Ret;
}
}
DGE_MATRICE_END