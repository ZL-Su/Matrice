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
#include "core/vector.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty> class _Quaternion {
	using _Myt = _Quaternion<_Ty>;
public:
	using value_type = _Ty;
	using vector_type = Vec3_<value_type>;

	_Quaternion() noexcept {
	}
	_Quaternion(value_type re, const vector_type& im) noexcept
		:m_real(re), m_imag(im) {
	}
	~_Quaternion() {
	}

	/**
	 *\brief accessor for the real part of the quaternion.
	 */
	MATRICE_GLOBAL_INL const value_type& real() const noexcept {
		return (m_real);
	}
	MATRICE_GLOBAL_INL value_type& real() noexcept {
		return (m_real);
	}

	/**
	 *\brief accessor for the imaginary part of the quaternion.
	 */
	MATRICE_GLOBAL_INL const vector_type& imag() const noexcept {
		return (m_imag);
	}
	MATRICE_GLOBAL_INL vector_type& imag() noexcept {
		return (m_imag);
	}

	/**
	 *\brief returns the norm of the quaternion.
	 */
	MATRICE_GLOBAL_INL value_type norm() const noexcept {
		return sqrt(sqr(m_real) + m_imag.dot(m_imag));
	}

	/**
	 *\brief returns the complex conjugate of the quaternion.
	 */
	MATRICE_GLOBAL_INL _Myt conj() const noexcept {
		return _Myt(m_real, 0 - m_imag);
	}

	/**
	 *\brief converts the quaternion to a rotation matrix.
	 */
	MATRICE_GLOBAL_INL Matrix_<value_type, 3> matrix() const noexcept {
		const auto s = 2*this->norm();
		Matrix_<value_type, 3> _Ret;
		_Ret[0][0] = 1 - s * (sqr(m_imag.y) + sqr(m_imag.z));
		_Ret[1][1] = 1 - s * (sqr(m_imag.x) + sqr(m_imag.z));
		_Ret[2][2] = 1 - s * (sqr(m_imag.x) + sqr(m_imag.y));
		_Ret[0][1] = s * (m_imag.x*m_imag.y - m_imag.z*m_real);
		return (_Ret);
	}

private:
	value_type  m_real = one<value_type>;
	vector_type m_imag;
};
_DETAIL_END
DGE_MATRICE_END