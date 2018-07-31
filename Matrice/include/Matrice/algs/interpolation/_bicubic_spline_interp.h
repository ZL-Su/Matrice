/**************************************************************************
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
**************************************************************************/
#pragma once

#include "_base.h"

MATRICE_ALGS_BEGIN
template<typename _Ty, size_t _Options = INTERP|BICUBIC|BSPLINE>
class BicubicSplineInterp : public InterpBase_<_Ty, BicubicSplineInterp<_Ty, _Options>>
{
	using base_t = InterpBase_<_Ty, BicubicSplineInterp<_Ty, _Options>>;
public:
	enum { Options = _Options };
	using typename base_t::value_t;
	using typename base_t::matrix_t;

	MATRICE_GLOBAL_FINL BicubicSplineInterp(const matrix_t& _Data)
		:base_t(_Data) {}

	MATRICE_HOST_INL void _Bspline_coeff(const matrix_t& _Data);

private:
	using base_t::m_coeff;
	const types::Matrix_<value_t, 4, 4> m_icoef{1, -3, 3, -1, 4, 0, -6, 3, 1, 3, 3, -3, 0, 0, 0, 1};
	const types::Matrix_<value_t, 4, 3> m_gcoef{-3, 6, -3, 0, -12, 9, 3, 6, -9, 0, 0, 3};
};
MATRICE_ALGS_END