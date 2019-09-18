/***********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include "_base.h"

MATRICE_ALGS_BEGIN
template<typename _Ty>
class _Bilinear_interpolation
	: public _Interpolation_base<_Bilinear_interpolation<_Ty>>
{
	using _Myt = _Bilinear_interpolation;
	using _Mybase = _Interpolation_base<_Myt>;
public:
	using typename _Mybase::matrix_type;
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;

	using _Mybase::_Interpolation_base;

	MATRICE_HOST_INL value_type operator()(const point_type& pos) const noexcept {
		return (*this)(pos.x, pos.y);
	}

	/**
	 * \brief computes the value at (x, y)
	 */
	MATRICE_HOST_INL value_type operator()(value_type x, value_type y) const noexcept {
		const auto x1 = max(floor<diff_t>(x), 0);
		const auto x2 = min(x1 + 1, _Mydata.cols() - 1);
		const auto y1 = max(floor<diff_t>(y), 0);
		const auto y2 = min(y1 + 1, _Mydata.rows() - 1);
		const auto f11 = _Mydata[y1][x1];
		const auto f12 = _Mydata[y2][x1];
		const auto f21 = _Mydata[y1][x2];
		const auto f22 = _Mydata[y2][x2];

		const auto dx1 = x - x1;
		const auto dx2 = x2 - x;
		return (y2 - y) * (dx2 * f11 + dx1 * f21) + (y - y1) * (dx2 * f12 + dx1 * f22);
	}

	MATRICE_HOST_INL auto grad(const point_type& pos) const noexcept {
		return (this->grad(pos.x, pos.y));
	}
	MATRICE_HOST_INL auto grad(value_type x, value_type y) const noexcept{
		const auto x1 = max(floor<diff_t>(x), 0);
		const auto x2 = min(x1 + 1, _Mydata.cols() - 1);
		const auto y1 = max(floor<diff_t>(y), 0);
		const auto y2 = min(y1 + 1, _Mydata.rows() - 1);
		const auto f11 = _Mydata[y1][x1];
		const auto f12 = _Mydata[y2][x1];
		const auto f21 = _Mydata[y1][x2];
		const auto f22 = _Mydata[y2][x2];

		const auto dfdx = (f21 - f11) * (y2 - y) + (f22 - f12) * (y - y1);
		const auto dfdy = (f12 - f11) * (x2 - x) + (f22 - f21) * (x - x1);
		const auto deno = (x2 - x1) * (y2 - y1);
		return std::make_tuple(dfdx / deno, dfdy / deno);
	}

	MATRICE_HOST_INL void _Coeff_impl() noexcept = delete;

private:
	using _Mybase::_Mydata;
};
MATRICE_ALGS_END