/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
3D Vision and Photo-Mechanics.
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
#include "util/_exception.h"
#include "math/_primitive_funcs.hpp"
#include <functional>

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Mty>
MATRICE_GLOBAL_INL auto _Union(const _Mty& _Left, const _Mty& _Right) {
#ifdef MATRICE_DEBUG
	DGELOM_CHECK(_Left.shape() == _Right.shape(), 
		"Shape inconsistent of the inputs in function _Union<_Mty>.");
#endif
	_Mty _Ret(_Left.shape());
	for (auto _Idx = 0; _Idx < _Left.size(); ++_Idx) {
		_Ret(_Idx) = max(_Left(_Idx), _Right(_Idx));
	}
	return forward<_Mty>(_Ret);
}
_DETAIL_END
DGE_MATRICE_END