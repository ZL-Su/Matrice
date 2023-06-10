/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2023, Zhilong(Dgelom) Su, all rights reserved.

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

template<class _Mty, typename _Cmp>
requires is_matrix_v<_Mty> || is_expression_v<_Mty>
MATRICE_HOST_INL auto _Argsort(const _Mty & _Arr, _Cmp && _Op) {
	Matrix_<size_t, _Arr.rows_at_compiletime, 1> _Ids(_Arr.rows());
	for (auto _Idx = 0; _Idx < _Ids.rows(); ++_Idx) {
		_Ids(_Idx) = _Idx;
	}
	MATRICE_USE_STD(stable_sort)(_Ids.begin(), _Ids.end(), _Op);
	return _Ids;
}
_DETAIL_END

template<typename... Args>
/**
 * \brief Sort an array by tracking the element index.
 * \param Args... Argpack of a `dgelom::Matrix` and a user defined comparator. 
          See `detail::_Argsort` for more details.
 * \returns Sorted indices in a linear `dgelom::Matrix`.
 */
MATRICE_HOST_FINL auto argsort(Args&&... _Args) noexcept {
	return detail::_Argsort(_Args...);
}
DGE_MATRICE_END
