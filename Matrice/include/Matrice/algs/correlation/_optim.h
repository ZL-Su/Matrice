/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
***********************************************************************/
#pragma once
#include "../../core/matrix.h"
#include "../../core/vector.h"
#include "../../core/solver.h"
#include "../../core/tensor.h"
#include "../../private/math/_linear.h"

MATRICE_ALGS_BEGIN namespace corr { _DETAIL_BEGIN

template<typename _Itptag> struct _Corr_border_size {};
template<> struct _Corr_border_size<_TAG bicspl_tag> {
	static constexpr auto lower = 1, upper = 2;
};
template<> struct _Corr_border_size<_TAG biqspl_tag> {
	static constexpr auto lower = 2, upper = 3;
};
template<> struct _Corr_border_size<_TAG bisspl_tag> {
	static constexpr auto lower = 3, upper = 4;
};

struct _Correlation_options {
	std::size_t _Stride =  7; //node spacing
	std::size_t _Radius = 10; //patch radius
	std::size_t _Maxits = 20; //maximum iterations

	template<typename _Ty>
	static constexpr _Ty _Mytol = _Ty(1.0e-6); //iteration tolerance

	/**
	 * \retrieve the range of patch centered on point _Pos
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL auto range(const _Ty& _Pos) {
		const auto x = floor<int>(_Pos.x), y = floor<int>(_Pos.y);
		return std::make_tuple(x - _Radius, x + _Radius, y - _Radius, y + _Radius);
	}
};

_DETAIL_END } MATRICE_ALGS_END