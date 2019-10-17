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
	 *\brief rotate a point with the quaternion.
	 *\param [point] a given 3d point.
	 */
	MATRICE_GLOBAL_INL vector_type rot(const vector_type& point) const noexcept {
		//_Myt p(0, point);
		auto prod = (*this)*_Myt(0, point)*this->conj();
		return (prod.imag());
	}

	/**
	 *\brief converts the quaternion to a rotation matrix.
	 */
	MATRICE_GLOBAL_INL Matrix_<value_type, 3> matrix() const noexcept {
		const auto s = 2*this->norm();

		Matrix_<value_type, 3> _Ret;
		const auto x = m_imag.x, y = m_imag.y, z = m_imag.z;
		_Ret[0][0] = 1 - s * (sqr(y) + sqr(z));
		_Ret[1][1] = 1 - s * (sqr(x) + sqr(z));
		_Ret[2][2] = 1 - s * (sqr(x) + sqr(y));
		_Ret[0][1] = s * (x*y - z * m_real);
		_Ret[0][2] = s * (x*z + y * m_real);
		_Ret[1][0] = s * (x*y + z * m_real);
		_Ret[1][2] = s * (z*y - x * m_real);
		_Ret[2][0] = s * (x*z - y * m_real);
		_Ret[2][1] = s * (y*z + x * m_real);

		return (_Ret);
	}

	MATRICE_GLOBAL_INL _Myt& operator+=(const value_type _real) noexcept {
		this->real() += _real;
		return (*this);
	}
	MATRICE_GLOBAL_INL _Myt& operator+=(const _Myt& _other) noexcept {
		this->real() += _other.real();
		this->imag() = m_imag + _other.imag();
		return (*this);
	}
	MATRICE_GLOBAL_INL _Myt& operator*=(const _Myt& _right) noexcept {
		const auto dot = m_imag.dot(_right.imag());
		const auto cro = m_imag.cross(_right.imag());
		m_imag = m_real * _right.imag() + _right.real()*m_imag + cro;
		m_real = m_real * _right.real() - dot;
		return (*this);
	}

	template<typename _Ty> friend
	MATRICE_GLOBAL_INL auto operator+(const _Quaternion<_Ty>& _left, const _Ty _real)noexcept {
		auto _Ret = _left; return (_Ret += _real);
	}
	template<typename _Ty> friend
	MATRICE_GLOBAL_INL auto operator+(const _Quaternion<_Ty>& _left, const _Quaternion<_Ty>& _right)noexcept {
		auto _Ret = _left; return (_Ret += _right);
	}
	template<typename _Ty> friend
	MATRICE_GLOBAL_INL auto operator*(const _Quaternion<_Ty>& _left, const _Quaternion<_Ty>& _right)noexcept {
		auto _Ret = _left; return (_Ret *= _right);
	}

private:
	value_type  m_real = one<value_type>;
	vector_type m_imag;
};
_DETAIL_END
DGE_MATRICE_END