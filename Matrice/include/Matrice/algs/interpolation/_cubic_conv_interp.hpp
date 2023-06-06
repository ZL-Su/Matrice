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
class _Cubic_conv_interpolation {
	using _Myt = _Cubic_conv_interpolation;
public:
	using value_type = _Ty;
	template<typename _Pixty>
	_Cubic_conv_interpolation(const Matrix<_Pixty>& img)
		:m_data(img.rows() + 1 + 2, img.cols() + 1 + 2, 0) {
		m_data.block(1, img.cols() + 1, 1, img.rows() + 1) = img;
		if constexpr (is_not_same_v<_Ty, _Pixty> && is_same_v<_Pixty, uint8_t>) {
			m_data = m_data / value_type(255);
		}
	}

	MATRICE_HOST_INL value_type operator()(value_type x, value_type y) const noexcept {
		const auto ix = floor<index_t>(x), iy = floor<index_t>(y);
		const auto dx = x - ix, dy = y - iy;
		const auto blk = m_data.block(ix, ix + 3, iy, iy + 3);

		value_type b[] = { 
			_Kernel(dx, blk[0]), _Kernel(dx, blk[1]),
			_Kernel(dx, blk[2]), _Kernel(dx, blk[3]) };

		return _Kernel(dy, b);
	}

	MATRICE_HOST_INL auto grad(value_type x, value_type y) const noexcept {

	}

private:
	MATRICE_HOST_INL value_type _Kernel(value_type dx, value_type* row) const noexcept {
		value_type v[] = {
				2 * row[1], row[2] - row[0],
				2 * row[0] - 5 * row[1] + 4 * row[2] - row[3],
				row[3] - row[0] + 3 * (row[1] - row[2])
		};
		return value_type(0.5) * (v[0] + v[1] * dx + v[2] * dx * dx + v[3] * dx * dx * dx);
	}
	Matrix<_Ty> m_data; //Default is zero-padding
};
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using bicubic_conv_interp = _Cubic_conv_interpolation<_Ty>;
MATRICE_ALGS_END