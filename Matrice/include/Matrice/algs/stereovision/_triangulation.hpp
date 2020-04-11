#pragma once
/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#include "core.hpp"
#include "../geometry/transform.h"

MATRICE_ALGS_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
class _Trig_base {
	using _Myt = _Trig_base;
public:
	using value_type = _Ty;
	using vector_type = Vec3_<value_type>;
	using matrix_2x4t = Matrix_<value_type, 2, 4>;
	using matrix_3x3t = Matrix_<value_type, 3, 3>;

	MATRICE_HOST_INL _Myt& set_extpars(initlist<value_type> ext) noexcept {
		decltype(auto) _It = ext.begin();
		const auto rx = *_It, ry = *(_It + 1), rz = *(_It + 2);
		const auto tx = *(_It + 3), ty = *(_It + 4), tz = *(_It + 5);

		rodrigues(vector_type{ rx, ry, rz }, m_rotmat);
		m_transv.x = tx, m_transv.y = ty, m_transv.z = tz;
		
		return (*this);
	}
	MATRICE_HOST_INL _Myt& set_intpars(nested_initlist<value_type> kk) noexcept {
		m_inpars.rview(0) = *kk.begin();
		m_inpars.rview(1) = *(kk.begin() + 1);
		return (*this);
	}

protected:
	matrix_2x4t m_inpars;
	matrix_3x3t m_rotmat;
	vector_type m_transv;
};
_DETAIL_END
MATRICE_ALGS_END