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
***********************************************************************/
#pragma once
#include "matrix.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
/// <summary>
/// \brief CLASS TEMPLATE, Smooth image type
/// </summary>
/// <typeparam name="_Ty">Floating point type</typeparam>
/// <typeparam name="_Ip">Interpolation operator type</typeparam>
template<typename _Ty, class _Ip> 
requires is_floating_point_v<_Ty>
class _Image : public Matrix_<_Ty, ::dynamic, ::dynamic> {
	using _Mybase = Matrix_<_Ty, ::dynamic, ::dynamic>;
	using _Myt = _Image;
	using _Myop = _Ip;
public:
	template<typename _Uy>
	requires is_not_same_v<_Ty, _Uy>
	_Image(Matrix_<_Uy, ::dynamic>&& _Src) {

	}

	/**
	 * \brief OP, eval image value at coordinates (x, y). 
	 */
	template<typename _Uy>
	decltype(auto) operator()(_Uy x, _Uy y) const {

	}
	template<typename _Uy>
	decltype(auto) operator()(_Uy x, _Uy y) {

	}

	/**
	 * \brief FUNCTION, serialize an '_Image<_Ty, _Ip>' object.
	 */
	template<typename _Ty, class _Ip>
	friend auto serialize(const _Image<_Ty, _Ip>& img) noexcept {

	}
private:
	_Myop _Myinst;
};

/// <summary>
/// \brief CLASS TEMPLATE, General image type
/// </summary>
/// <typeparam name="_Ty">Scalar value type</typeparam>
template<typename _Ty>
class _Image : public Matrix_<_Ty, ::dynamic, ::dynamic> {

};
_DETAIL_END
DGE_MATRICE_END