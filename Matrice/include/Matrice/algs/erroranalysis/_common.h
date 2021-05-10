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
**********************************************************************/
#pragma once

#include "core/matrix.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
class _Error_analysis {
	using value_type = _Ty;
	using array_type = Matrix<_Ty>;
public:
	MATRICE_HOST_INL 
	static auto mean_bias(const array_type& _Vals, const value_type _Ref) {
		return (_Vals - _Ref).sum() / _Vals.size();
	}

};
_DETAIL_END
DGE_MATRICE_END