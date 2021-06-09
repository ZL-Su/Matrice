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
}
DGE_MATRICE_END