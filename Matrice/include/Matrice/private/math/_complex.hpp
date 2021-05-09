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
#include "core/matrix.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty, size_t _M, size_t _N> class _Complex{};

/// <summary>
/// \brief Partial specialization
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty, size_t _M, size_t _N>
requires is_scalar_v<_Ty> 
class _Complex<_Ty, _M, _N> 
{
public:
	using value_type = _Ty;
	using real_type = Matrix_<value_type, _M, _N>;
	using imag_type = real_type;
	using element_type = _Complex<value_type, 1, 1>;

	/**
	 * \brief Get complex element in a linear manner.
	 */
	MATRICE_HOST_INL element_type operator[](index_t i) const noexcept {
		return element_type{ _Myreal(i), _Myimag(i) };
	}
	/**
	 * \brief Get complex element in a grid manner.
	 */
	MATRICE_HOST_INL element_type operator()(index_t i, index_t j) const noexcept {
		return element_type{ _Myreal[i][j], _Myimag[i][j] };
	}

private:
	real_type _Myreal;
	imag_type _Myimag;
};

/// <summary>
/// \brief Specialization to primary complex number: z = _Myreal + i*_Myimag.
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty>
requires is_scalar_v<_Ty>
class _Complex<_Ty, 1, 1> 
{
	using _Myt = _Complex<_Ty, 1, 1>;
public:
	using value_type = _Ty;
	using real_type = value_type;
	using imag_type = real_type;

	_Complex() = default;
	_Complex(value_type&& _Re, value_type&& _Im) noexcept
		:_Myreal{ _Re }, _Myimag{ _Im } {
	}
	_Complex(const value_type& _Re, const value_type& _Im) noexcept
		:_Myreal{ _Re }, _Myimag{ _Im } {
	}

	/**
	 * \brief Get real part of the complex number.
	 */
	MATRICE_GLOBAL_FINL decltype(auto) real() const noexcept {
		return (_Myreal);
	}
	MATRICE_GLOBAL_FINL decltype(auto) real() noexcept {
		return (_Myreal);
	}

	/**
	 * \brief Get imaginary part of the complex number.
	 */
	MATRICE_GLOBAL_FINL decltype(auto) imag() const noexcept {
		return (_Myimag);
	}
	MATRICE_GLOBAL_FINL decltype(auto) imag() noexcept {
		return (_Myimag);
	}

	/**
	 * \brief Get binded real and imaginary parts.
	 * \return auto [re, im] = complex;
	 */
	MATRICE_GLOBAL_FINL operator auto() const noexcept {
		return MATRICE_STD(tuple) { _Myreal, _Myimag };
	}

	/**
	 * \brief Get this complex number. 
	   Compatible with the leading _Complex<_Ty, _M, _N>.
	 */
	MATRICE_GLOBAL_FINL decltype(auto) operator[](auto) const noexcept {
		return (*this);
	}
	MATRICE_GLOBAL_FINL decltype(auto) operator()(auto...) const noexcept {
		return (*this);
	}

	/**
	 * \brief Get the conjugate number.
	 */
	MATRICE_GLOBAL_FINL _Myt conj() const noexcept {
		return _Myt{_Myreal, -_Myimag };
	}

	/**
	 * \brief Complex number multiplication.
	 */
	MATRICE_GLOBAL_FINL _Myt operator*(const _Myt& _Other) const noexcept {
		const auto _Real = real() * _Other.real() - imag() * _Other.imag();
		const auto _Imag = real() * _Other.imag() + imag() * _Other.real();
		return _Myt{_Real, _Imag};
	}

	/**
	 * \brief Complex number multiplies a real number.
	 */
	MATRICE_GLOBAL_FINL _Myt operator*(value_type _Ampl) const noexcept {
		return _Myt{ _Ampl*real(), _Ampl*imag() };
	}
	MATRICE_GLOBAL_FINL
	friend _Myt operator*(_Ty _Left, const _Myt& _Right) noexcept {
		return _Right.operator*(_Left);
	}

public:
	real_type _Myreal{0};
	imag_type _Myimag{0};
};

_DETAIL_END

/// <summary>
/// \brief ALIAS TEMPLATE, generic type for complex number or structure.
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty, size_t _M, size_t _N = _M>
requires is_scalar_v<_Ty>
using complex_t = detail::_Complex<_Ty, _M, _N>;

DGE_MATRICE_END