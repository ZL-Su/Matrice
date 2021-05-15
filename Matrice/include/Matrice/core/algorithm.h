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
#include <algorithm>

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
/**
 * INTERNAL FUNC, get the index of the maximum value in range [_First, _Last).
 */
template<typename _FwdIt>
MATRICE_DEVICE_INL size_t _Argmax(_FwdIt _First, _FwdIt _Last) noexcept {
	auto _Iter = _First;
	auto _Max = _Iter;
	for (; _Iter != _Last; ++_Iter) {
		if (*_Max < *_Iter) {
			_Max = _Iter;
		}
	}
	return (_Max - _First);
}
_DETAIL_END
DGE_MATRICE_END
